  !========================================================================
  !             subroutine new_fem_solver
  !========================================================================
  subroutine new_fem_solver( this, prefix, mat_type )
    use precision,only:wp,doubleprec
    use global_parameters, only: mype,mtheta,mpsi,numberpe
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscmat.h" 
#include "finclude/petscvec.h" 
#include "finclude/petscviewer.h"
    type(solver_data) :: this
    character*(*) prefix
    integer mat_type
    ! locals
    real(wp) :: alpha,beta
    Mat Amat, Mmat
    Mat Bmat
    integer ierr,neq,mindex,mring, i, j, ir, it, ib, ib2, id
    PetscScalar a
    PetscViewer viewer
    integer, dimension(:), allocatable :: bdy
    Vec              vec_yy

    interface
       subroutine form_tangent(Amat, Mmat, iemax, iet1, iet2, iet3, npts, points, &
            alpha, beta, ierr)
         use precision,only:wp,doubleprec
         implicit none
         integer iemax,nptsd,npts,ierr
         integer,dimension(:) :: iet1,iet2,iet3
         real(wp) points(2,npts)
         real(wp) alpha, beta
         Mat Amat, Mmat
       end subroutine form_tangent
    end interface

    mring=2
    mindex=64*mring+1 ! estimate for matrix size
    neq = sum( mtheta )

!cj    call solver_create_mat( this, Amat, 10, neq, ierr )
!cj    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
!cj    call MatSetOption( Amat, MAT_IGNORE_OFF_PROC_ENTRIES, ierr )
    call MatCreate( partd_comm, Amat, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)    
    call MatSetSizes(Amat, PETSC_DECIDE,PETSC_DECIDE, neq, neq, ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatSetFromOptions( Amat, ierr ) ! need ti separate for SuperLU
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex/4+1,PETSC_NULL_INTEGER,ierr)
    call MatSeqAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    call MatCreate( partd_comm, Mmat, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)    
    call MatSetSizes(Mmat, PETSC_DECIDE,PETSC_DECIDE, neq, neq, ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatSetFromOptions( Mmat, ierr ) ! need ti separate for SuperLU
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatMPIAIJSetPreallocation(Mmat,mindex,PETSC_NULL_INTEGER,mindex/4+1,PETSC_NULL_INTEGER,ierr)
    call MatSeqAIJSetPreallocation(Mmat,mindex,PETSC_NULL_INTEGER,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    alpha = 1.
    beta  = 1.
    call form_tangent(Amat,Mmat,this%fem%iemax,this%fem%iet1,this%fem%iet2,&
         this%fem%iet3,neq,this%fem%points,alpha,beta,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
!cj  do i=1,neq
!cj        call MatSetValue(Amat,1,i,1,i,alpha,ADD_VALUES,ierr)
!cj  enddo

    call MatAssemblyBegin( Amat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyEnd( Amat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyBegin( Mmat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyEnd( Mmat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
!cj    if(mype==0) write(*,*) " ---> before bc fix", ierr
!cj       call PetscViewerBinaryOpen(partd_comm, 'A0.m',FILE_MODE_WRITE,viewer,ierr)
!cj       call MatView( Amat, viewer, ierr )
!cj       call PetscViewerDestroy(viewer,ierr)
    if(mype==0) write(*,*) "	finite elements matrix created"

    call new_gk_matid( this%mat, mat_type )

!cj imposing dirichlet bc
    if( mat_type .eq. 1 ) ib = mtheta(0)
    if( mat_type .eq. 2 ) ib = mtheta(0) + mtheta(mpsi)
    allocate(bdy(0:ib-1))
    this%mat%numbdy = mtheta(0) + mtheta(mpsi)
    allocate(this%mat%bdy(0:ib-1))
    it = 0
    ir = 0
    id = 0
    ib2 = 0
    allocate(this%zero01(neq),this%zero10(neq))
    this%zero01=0.
    this%zero10=1.
    do i=0,mpsi
       do j=0,mtheta(i)
          it = it + 1 !global index
          if( j == mtheta(i) ) then
             ir = ir + 1   !repeated index
          else
#ifdef TODO
             if( mat_type .eq. 1 .and. i==0 ) then
                ! use 0-based bdy index to set Dirichlet bc for lhs matrix
                bdy(ib2)=id
                this%mat%bdy(ib2)=id
                ib2=ib2+1  !count bdy index and assign global index

                ! use 1-based bdy index to set Dirichlet bc for rhs vector
                this%zero01(1+id) = 1.
                this%zero10(1+id) = 0.
!cj    if(mype .eq. 0 ) write(112,112) 1+id, this%zero01(1+id), this%zero10(1+id)
!cj112 format("finite element solver testsolve bdy=", 1(2x,i7), 2(1x,d12.6))
             endif
             if( mat_type .eq. 2 .and. (i==0 .or. i==mpsi) ) then
                ! use 0-based bdy index to set Dirichlet bc for lhs matrix
                bdy(ib2)=id
                this%mat%bdy(ib2)=id
                ib2=ib2+1  !count bdy index and assign global index

                ! use 1-based bdy index to set Dirichlet bc for rhs vector
                this%zero01(1+id) = 1.
                this%zero10(1+id) = 0.
!cj    if(mype .eq. 0 ) write(112,112) 1+id, this%zero01(1+id), this%zero10(1+id)
!cj112 format("finite element solver testsolve bdy=", 1(2x,i7), 2(1x,d12.6))
             endif
#else
              if ( this%mat%gk_matid(it) .lt. 0 )  then
                 ! use 0-based bdy index to set Dirichlet bc for lhs matrix
                 bdy(ib2)=id
                 this%mat%bdy(ib2)=id
                 ib2=ib2+1  !count bdy index and assign global index
 
                 ! use 1-based bdy index to set Dirichlet bc for rhs vector
                 this%zero01(1+id) = 1.
                 this%zero10(1+id) = 0.
              endif
#endif
             id = id + 1   !conunt global index
          end if
       end do
    end do
    if(mype==0) write(*,*) "	finite elements bdy ib ib2=", ib, ib2
    if(ib .ne. ib2) call MPI_ABORT(mpi_comm_world,1,ierr)
!cj    if(mype==0) write(113,111) ir, ib, mtheta(0), mtheta(mpsi)
!cj111 format("finite element solver bdy index =", 4(2x,i7))
!cj    do i = 0, ib-1
!cj    if(mype==0) write(113,113) i, bdy(i)
!cj    enddo
!cj113 format("bdy index =", 2(2x,i7))

!   Set Amat=-1 on bdy vertex for applying Dirichlet bc
!   using zero-based bdy index
    a = -1.
    call MatZeroRows(Amat,ib,bdy,a,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
    if(mype==0) write(*,*) "	finite elements matrix bdy fix"
!cj MatZeroRows(Mat, numRows,const PetscInt rows[], diag,Vec x,Vec b)
!cjISCreateGeneral(poloidalComm[tRank],
!cj                NumBoundaryVertices, bdyRowIndex,
!cj                PETSC_COPY_VALUES, isBdy);
!cjMatZeroRows(userSolver.lhsMatrix, map.NumBoundaryVertices,
!cj            bdyRowIndex, a /*Sone*/

!cj       call PetscViewerBinaryOpen(partd_comm, 'A1.m',FILE_MODE_WRITE,viewer,ierr)
!cj       call MatView( Amat, viewer, ierr )
!cj       call PetscViewerDestroy(viewer,ierr)

    call solver_init( this, prefix, Amat, neq,ierr)
    if(mype==0) write(0,*) '	finite elements solver created'
    
    ! finalize stiff matrix
    ! ... done in solver_init

    ! finalize mass matrix
    call MatGetSubMatrix(Mmat,this%petsc_gts_is,this%petsc_gts_is,MAT_INITIAL_MATRIX,Bmat,ierr)
    Mmat=Bmat
    this%mat%mdata = Mmat

    ! create zeroVec01 for applying bc
    call VecDuplicate(this%bVec, this%zeroVec01, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecCreateSeqWithArray(MPI_COMM_SELF,neq,this%zero01,vec_yy,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterBegin(this%scat,vec_yy, this%zeroVec01, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(this%scat,vec_yy, this%zeroVec01, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

    ! create zeroVec10 for applying bc
    call VecDuplicate(this%bVec, this%zeroVec10, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecCreateSeqWithArray(MPI_COMM_SELF,neq,this%zero10,vec_yy,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterBegin(this%scat,vec_yy, this%zeroVec10, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(this%scat,vec_yy, this%zeroVec10, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    call VecDestroy( vec_yy, ierr )

!cj    ! scale smoother
!cj    call MatGetDiagonal(Amat,this%bVec,ierr)
!cj    call VecReciprocal(this%bVec,ierr)
!cj    ! ajust weight: weight for myself is (1+beta)
!cj    ! beta=-.5d0
!cj    beta=-.375d0
!cj    call VecScale(this%bVec,beta,ierr)
!cj    call MatDiagonalScale(Amat,this%bVec,PETSC_NULL,ierr)
!cj    beta = 1.d0
!cj    call MatShift(Amat,beta,ierr)

!cj! debug to check matrix
!cj    if( mype / nproc_partd .eq. -10 ) then
!cj       call PetscViewerASCIIOpen(partd_comm, 'A.m',viewer,ierr)
!cj       call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr)
!cj       call MatView( Mmat, viewer, ierr )
!cj       call PetscViewerDestroy( viewer, ierr )
!cj    endif

    ! clear points here - not needed any more
    deallocate(bdy)
    deallocate(this%fem%points)
    nullify(this%fem%points)
    deallocate(this%fem%iet1)
    nullify(this%fem%iet1)
    deallocate(this%fem%iet2)
    nullify(this%fem%iet2)
    deallocate(this%fem%iet3)
    nullify(this%fem%iet3)
    if ( this%fem%order .gt. 1 ) then
       deallocate(this%fem%points_o2)
       nullify(this%fem%points_o2)
       deallocate(this%fem%iet4)
       nullify(this%fem%iet4)
       deallocate(this%fem%iet5)
       nullify(this%fem%iet5)
       deallocate(this%fem%iet6)
       nullify(this%fem%iet6)
    endif

  end subroutine new_fem_solver


!========================================================================
!         testsolve
!========================================================================
  subroutine testsolve( solver )
    !     solve linear system.  
    !     Solver (KSP) and matrix are in linear_solve.h
    !     
    !     input:
    !       nn: global number of equations in system
    !       bb_gk: rhs[nn]
    !       xx_gk: solution[nn]
    !     
    !     output:
    !       bb_gk: RHS[nn]
    !     ierr: error code
    !     
    !use global_parameters
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mtheta,mpsi,mgrid,zetamin,deltaz,deltat,R_0
    use field_array,only: psi_grid,qtinv
    use esi_interface
    implicit none
    type(solver_data) :: solver
    integer :: nn
    integer :: ierr
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq
    real(doubleprec),allocatable,dimension(:) :: xx, rhs 
!cj    real(doubleprec),allocatable,dimension(:) :: zero01, zero10 
    integer count
    real(wp) zeta, pdum, dt, tdum, tmp
    real(wp), dimension(2) :: aaa,the,rj,zj,Bj
    real(wp) M, L, alfa, PI
    data count/0/
    count = count + 1

    write(0, *) "	finite element test solve start"

!cj the following is a combianation of solve and tri_mesh
    neq = sum(mtheta)
    allocate(xx(neq),rhs(neq))
    xx=0.
    ivsl = 0
    ivgts = 0
    zeta = zetamin + deltaz ! zeta for my plane   
    do ir=0,mpsi 
       itmax=mtheta(ir)
       pdum=psi_grid(ir)
       dt = deltat(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             tdum = real(it)*dt + zeta*qtinv(ir)
             aaa(1)=pdum
             the(1)=tdum
             call esigetrzb(rj,zj,Bj,aaa,the,1)
!cj u = x^2 + y^2, del u = 4
!cj             rhs(ivsl) = 4.
             M=1.
             L=1.
             alfa=0.
             PI=3.14159265359
             xx(ivsl) = sin(PI*M*(rj(1)-R_0)) * cos(L*zj(1)+alfa)
             rhs(ivsl) = (-(PI*M)*(PI*M) - L*L) * xx(ivsl) 
             !cj if( ir==0 .or. ir==mpsi ) then
             !cj must be set to negative xx, since matrix was set to negative
!cj u = x^2 + y^2, negative since a=-1 set as Amat diagonal
!cj                rhs(ivsl) = - xx(ivsl)
!cj u = x^2 + y^2, del u = 4
!cj             xx(ivsl) = -( (rj(1)-R_0)*(rj(1)-R_0) + zj(1)*zj(1) )
             ! assign the boundary vertex
             if( solver%zero01(ivsl) .gt. 0 ) then
                xx(ivsl) = - sin(PI*M*(rj(1)-R_0)) * cos(L*zj(1)+alfa)
             endif
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    
    ! solve GK Poisson equation
    itestsolve=333
    call petsc_solve( solver, neq, rhs, xx, ierr)

!cj the following is a combianation of solve and tri_mesh, check |u-xx|
    ivsl = 0
    ivgts = 0
    zeta = zetamin + deltaz ! zeta for my plane   
    do ir=0,mpsi 
       itmax=mtheta(ir)
       pdum=psi_grid(ir)
       dt = deltat(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1

             tdum = real(it)*dt + zeta*qtinv(ir)
             aaa(1)=pdum
             the(1)=tdum
             call esigetrzb(rj,zj,Bj,aaa,the,1)

!cj u = x^2 + y^2
!cj             tmp = ( (rj(1)-R_0)*(rj(1)-R_0) + zj(1)*zj(1) )
!cj    if(mype .eq. 0 .and. ir.ne.0 .and. ir.ne.mpsi ) write(111,111) ivsl, xx(ivsl), tmp
             tmp = sin(PI*M*(rj(1)-R_0)) * cos(L*zj(1)+alfa)
             ! print the inner vertex
             !if( mype.eq.0 .and. solver%zero10(ivsl) .gt.0 ) write(111,111) ivsl, rj(1)-R_0, zj(1), xx(ivsl), tmp
             if( mype.eq.0 .and. solver%zero10(ivsl) .gt.0 ) write(111,111) ivsl, xx(ivsl), tmp
111 format("finite element solver testsolve solution=", 1(2x,i7), 2(1x,d12.6))
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    call flush(111)
    write(0, *) "	finite element test solve end"
    
    deallocate(rhs,xx) 

    ! debug
    if( count .eq. 1 .and. mype==-10 ) then
       call plot_array( neq, xx, 'x.m' )
       call plot_array( neq, rhs, 'b.m' )
       stop 'plot'
    endif

  end subroutine testsolve


!========================================================================
!         petsc_solve
!========================================================================
  subroutine petsc_solve( solver, nn, bb, xx, ierr )
    !     solve linear system.  
    !     Solver (KSP) and matrix are in linear_solve.h
    !     
    !     input:
    !       nn: global number of equations in system
    !       bb: rhs[nn]
    !       xx: solution[nn]
    !     
    !     output:
    !       bb: RHS[nn]
    !     ierr: error code
    !     
    !use global_parameters
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mgrid ! debug
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
    integer :: nn,ierr
    real(doubleprec),intent(in) :: bb(nn)
    real(doubleprec),intent(inout) :: xx(nn)
    type(solver_data) :: solver
    ! locals
    Vec              vec_yy, rhsVec, tmpVec01, tmpVec10
    integer          iplot,i,j
    PetscScalar t1
!    PetscInt,allocatable,dimension(:) :: idx
    PetscScalar,allocatable,dimension(:) :: yy, zz
!    PetscInt,allocatable,dimension(:) :: idx
    PetscViewer viewer
    integer :: m,n,p
    PetscInt :: its
    PetscReal :: rnorm

    if(solver%mat%data .eq. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    !     copy real vector into double precision
    allocate(yy(nn),zz(nn))
!    allocate(idx(nn))

    !     scatter rhs into PETSc vector
    yy = bb
    call VecCreateSeqWithArray(MPI_COMM_SELF,nn,yy,vec_yy,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    !     call VecSet( solver%bVec, 1.e30, ierr ) ! debug
    call VecScatterBegin(solver%scat,vec_yy, solver%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(solver%scat,vec_yy, solver%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

    if(itestsolve==333) then 
    !cj M * b =>> rhs for solve Stiff xx = rhs
    call VecDuplicate(solver%bVec, rhsVec, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    call VecCopy(solver%bVec, rhsVec, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    call MatMult(solver%mat%mdata, rhsVec, solver%bVec, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    endif

    if(itestsolve==333) then 
!    do i = 1, nn !!solver%mat%numbdy
!    if(mype==0) write(113,113) i, solver%zero10(i)
!!cj    if(mype==0) write(113,113) i, solver%zero01(i)
!    enddo
!113 format("bdy index =", 1(2x,i7), 1(1x, d12.6))


!   apply dirichlet bc
!   bVec = bVec * zero10 + xVec * zero01
    yy=xx
    call VecCreateSeqWithArray(MPI_COMM_SELF,nn,yy,vec_yy,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterBegin(solver%scat,vec_yy, solver%xVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(solver%scat,vec_yy, solver%xVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

!   take the part of bVec corresponding to internal vertex
    call VecDuplicate(solver%bVec, tmpVec10, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    call VecPointwiseMult(tmpVec10, solver%bVec, solver%zeroVec10, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

!   take the part of xVec corresponding to boundary vertex
    call VecDuplicate(solver%bVec, tmpVec01, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    call VecPointwiseMult(tmpVec01, solver%xVec, solver%zeroVec01, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

    call VecAXPBYPCZ(solver%bVec, 1., 1., 0., tmpVec10, tmpVec01, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    endif

!    if(itestsolve==444) then 
!    call VecSet(solver%xVec, 1., ierr)
!    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
!    call MatMult(solver%mat%data, solver%xVec, solver%bVec, ierr)
!    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
!    call VecSet(solver%xVec, 0., ierr)
!    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
!    endif
    
    !     Solve      
    call KSPSolve( solver%ksp, solver%bVec, solver%xVec, ierr )
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

!    if(itestsolve==444) then 
!    call KSPMonitor(solver%ksp, its, rnorm, ierr)
!    call KSPView(solver%ksp,PETSC_VIEWER_STDOUT_WORLD, ierr)
!    write(0,111) its, rnorm
!111 format("finite element testsolve its rnorm=", 1(2x,i7), 1(1x,d12.6))
!    if(mype==0) write(0,*) 'finite element solver testsolve xVec'
!    call VecView( solver%xVec, PETSC_VIEWER_STDOUT_WORLD, ierr )
!    endif

    !     scatter solution to gts array (yy)
    t1 = 0.0
    call VecSet( vec_yy, t1, ierr )
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterBegin(solver%scat,solver%xVec,vec_yy,INSERT_VALUES,SCATTER_REVERSE,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(solver%scat,solver%xVec, vec_yy, INSERT_VALUES,SCATTER_REVERSE,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr) 
    
    ! add all to all -- not super efficient
    call MPI_ALLREDUCE(yy,zz,nn,MPI_DOUBLE_PRECISION,MPI_SUM,partd_comm,ierr)
    !     copy in
    xx = zz;
    ! clean up
    call VecDestroy( vec_yy, ierr )
    call VecDestroy( rhsVec, ierr )
!    call VecDestroy( zeroVec01, ierr )
!    call VecDestroy( zeroVec10, ierr)
    call VecDestroy( tmpVec01, ierr )
    call VecDestroy( tmpVec10, ierr )
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)      
    deallocate(yy,zz)
!    deallocate(idx)
  end subroutine petsc_solve


!========================================================================
!         set up gk_matid for boundary condition
!========================================================================
  subroutine new_gk_matid( this, mat_type )
    use precision,only:wp,doubleprec
    use global_parameters,only: mype,mtheta,mpsi,ZF_SEPARATION,ifelectron
    use particle_decomp,only:myrank_partd,nproc_partd
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
    type(phit_matrix_data) :: this
    integer mat_type
    ! locals
    integer i, j, ierr, it, id

    ! form map from gts node numbering (includes 0 & 2pi) to matrix numbering
    it = sum(mtheta) + mpsi + 2
    allocate(this%gk_matid(it))
    this%gk_matid = -1
    this%gk_matid(it) = -2
    
    id = 0
    it = 0
    do i=0,mpsi
       do j=0,mtheta(i)
          it = it + 1
          if( j == mtheta(i) ) then
             this%gk_matid(it) = id - mtheta(i) + 1 - 1
          else
             id = id + 1
             this%gk_matid(it) = id - 1
          end if

          ! BCs
          if(mat_type==1) then      ! for Poisson solver
!cj             if(ZF_SEPARATION==1) then
                if( i==0 ) this%gk_matid(it) = -1         ! Dirichlet boundary condiction at inner boundary: phi(i=0)=0
!cj             else ! ZF_SEPARATION==0                
                if(ifelectron==0)  then ! adiabatic electrons
                   ! NOTE: any combination of following BCs produces a spike at the boundary for ad-electron case
                   ! if( i==0 ) this%gk_matid(it) = -1
                   ! if( i==mpsi ) this%gk_matid(it) = -1  ! Dirichlet boundary condiction at outer boundary: phi(i=mpsi)=0
                else
                   !   if( i==0 ) this%gk_matid(it) = -1
                   !   if( i==mpsi ) this%gk_matid(it) = -1  ! Dirichlet boundary condiction at outer boundary: phi(i=mpsi)=0
                endif
!cj             endif
          else                      ! for dphi/dt solver
             if( i==0 ) this%gk_matid(it) = -1
             if( i==mpsi ) this%gk_matid(it) = -1
          endif
          if(mype .eq. 0 ) write(112,112) it, this%gk_matid(it)
112 format("new_gk_matid ", 2(2x,i7))
       end do
    end do
    if(id .ne. sum(mtheta)) call MPI_ABORT(mpi_comm_world,1,ierr)
    if(this%gk_matid(sum(mtheta) + mpsi + 2) .ne. -2) call MPI_ABORT(mpi_comm_world,1,ierr)
  end subroutine new_gk_matid

end module solver_data_class


!========================================================================
!                         subroutine form_tangent
!========================================================================
subroutine form_tangent( Amat, Mmat, iemax, iet1, iet2, iet3, npts, points, &
     alpha, beta, ierr )
  !     
  !     input:
  !     iemax: number of elements
  !     iet1,2,3[iemax]: nodes of triangles of mesh
  !     npts: number of mesh points
  !     points[2,npts]: coordinates of points
  !     
  !     output:
  !     Amat: matrix to add tangent to
  !     ierr: error code
  !
  use precision,only:wp,doubleprec
  use particle_decomp,only:myrank_partd
  use global_parameters, only:mype
  implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
  integer iemax,npts,ierr
  integer,dimension(:) :: iet1,iet2,iet3
  real(wp) :: points(2,npts)
  real(wp) alpha, beta
  Mat  Amat, Mmat
  real(doubleprec)  xlt(2,3),sarr(3,3),sarr2(3*3),ul(3,4),eps ! therm2d uses *8
  real(doubleprec)           marr(3,3),marr2(3*3)             ! mass matrix
  real(doubleprec) c1, c2
  integer ie,i,j,k,kk,nz,ind2(3),ind3(3),ind(3),nel,nst

  ! Lagrange finite elements
  integer iw, is, il
  integer fe_order, fe_node, fe_ns, fe_nw
  real(doubleprec),dimension(:),allocatable :: fe_weight
  real(doubleprec),dimension(:,:),allocatable :: fe_NN
  real(doubleprec),dimension(:,:,:),allocatable :: fe_dN1, fe_Eij
  real(doubleprec) fe_mass(3,3),fe_stiff(3,3) ! mass and stiffness matrix

  c1 = alpha
  c2 = beta
  nel = 3                   ! points per element
  nst = 3                   ! size of elem matrix
  ul = 0.                   ! displacement, velocit and acceloration (not used)
  eps = 1.e-6

  call MatZeroEntries(Amat,ierr)

  call lagrange_element_define(fe_order, &
                               fe_node, &
                               fe_ns, &
                               fe_nw )
  allocate(fe_weight(fe_nw))
  allocate(fe_NN(fe_ns,fe_nw))
  allocate(fe_dN1(fe_node,fe_ns,fe_nw))
  allocate(fe_Eij(fe_ns,fe_ns,fe_nw))
  call lagrange_element_construct(fe_order, &
                                   fe_ns, &
                                   fe_nw, &
                                   fe_weight, &
                                   fe_NN, &
                                   fe_dN1, &
                                   fe_Eij)
!  do iw=1,fe_nw
!  if(mype .eq. 0) write(0,*) 'lagrange_element_construct: weight=', iw, fe_weight(iw)
!     do is=1,fe_ns
!  if(mype .eq. 0) write(0,*) 'lagrange_element_construct: NN=', is, fe_NN(is,iw)
!        do il=1,fe_node
!  if(mype .eq. 0) &
!  write(0,*) 'lagrange_element_construct: dN1=', &
!  iw, is, il, fe_dN1(il, is, iw)
!        enddo
!     do il=1,fe_ns
!  if(mype .eq. 0) &
!  write(0,*) 'lagrange_element_construct: Eij=', &
!  iw, is, il, fe_Eij(il, is, iw)
!     enddo
!     enddo
!  enddo

  if(mype .eq. 0) write(0,*) 'form_tangent: iemax=', iemax
  do ie=1,iemax
     ind(1)=iet1(ie)
     ind(2)=iet2(ie)
     ind(3)=iet3(ie)
     if(ind(1) .eq. 0) stop 'ind(1) == 0'
    
     xlt(1,1) = points(1,ind(1))
     xlt(1,2) = points(1,ind(2))
     xlt(1,3) = points(1,ind(3))
     xlt(2,1) = points(2,ind(1))
     xlt(2,2) = points(2,ind(2))
     xlt(2,3) = points(2,ind(3))

        sarr = 0.
        marr = 0.
        call therm2d(c1,c2,ul,xlt,sarr,marr,nel,nst,2,1)
        call lagrange_element_matrix(fe_ns,  &
                                     fe_nw,  &
                                     xlt(1,:),  &
                                     xlt(2,:),  &
                                     fe_weight,  &
                                     fe_NN, &
                                     fe_dN1, &
                                     fe_Eij, &
                                     fe_mass, &
                                     fe_stiff)

        nz = 0
        do i=1,nel
           k = ind(i)
           nz = nz + 1
           ind3(nz) = i
           ind2(nz) = k - 1 ! zero based
        enddo
!        if(mype==0) write(114,114) ie, nz, ind(1), ind(2), ind(3), &
!         fe_mass(1,1), fe_mass(1,2), fe_mass(1,3), &
!         fe_mass(2,1), fe_mass(2,2), fe_mass(2,3), &
!         fe_mass(3,1), fe_mass(3,2), fe_mass(3,3)
!114 format("matsetvalue index =", 5(2x,i7), 9(2x,d12.6))

        kk = 0
        do i=1,nel
           k = ind(i)
           do j=1,nz
             kk = kk + 1
             sarr2(kk) = sarr(i,ind3(j))
             marr2(kk) = marr(i,ind3(j))
!            sarr2(kk) = fe_stiff(i,ind3(j))
!            marr2(kk) = fe_mass(i,ind3(j))
           enddo
        enddo
        !
        call MatSetValues(Amat,nz,ind2,nz,ind2,sarr2,ADD_VALUES,ierr)
        call MatSetValues(Mmat,nz,ind2,nz,ind2,marr2,ADD_VALUES,ierr)
!        if(mype .eq. 0) write(0,*) 'form_tangent: ', ie,') nz=',nz,', ids=',ind2,' data=',arr2
  enddo

  deallocate(fe_weight)
  deallocate(fe_NN)
  deallocate(fe_dN1)
  deallocate(fe_Eij)
end subroutine form_tangent
