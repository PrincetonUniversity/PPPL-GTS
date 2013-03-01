!========================================================================
!                         module phit_matrix_type
!                           - shell to construct Phi tilda
!========================================================================
module phit_matrix_type
  implicit none
#include "finclude/petscdef.h"
#include "finclude/petscmatdef.h"
  type phit_matrix_data
     Mat :: data ! not allcatted here
     integer,dimension(:),pointer :: gk_matid
  end type phit_matrix_data
end module phit_matrix_type
!========================================================================
!                         module phit_matrix_class
!========================================================================
module phit_matrix_class
  use phit_matrix_type
  implicit none
  public :: add_gk_value,new_phit_matrix
contains
  !========================================================================
  !                       subroutine new_phit_matrix  
  !========================================================================
  subroutine new_phit_matrix( this, mat, mindex, mat_type )
    use precision,only:wp,doubleprec
    use global_parameters,only: mype,mtheta,mpsi,ZF_SEPARATION,ifelectron
    use particle_decomp,only:myrank_partd,nproc_partd
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
    type(phit_matrix_data) :: this
    Mat mat
    !integer,dimension(0:) :: mtheta
    integer mindex, mat_type
    ! locals
    integer i, j, ierr, it, id
    !real(doubleprec) val

    this%data = mat

    !if(it .ne. sum(mtheta)) call MPI_ABORT(mpi_comm_world,1,ierr)
    !end do
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
             if(ZF_SEPARATION==1) then
                if( i==0 ) this%gk_matid(it) = -1         ! Dirichlet boundary condiction at inner boundary: phi(i=0)=0
             else ! ZF_SEPARATION==0                
                if(ifelectron==0)  then ! adiabatic electrons
                   ! NOTE: any combination of following BCs produces a spike at the boundary for ad-electron case
                   ! if( i==0 ) this%gk_matid(it) = -1
                   ! if( i==mpsi ) this%gk_matid(it) = -1  ! Dirichlet boundary condiction at outer boundary: phi(i=mpsi)=0
                else
                   !   if( i==0 ) this%gk_matid(it) = -1
                   !   if( i==mpsi ) this%gk_matid(it) = -1  ! Dirichlet boundary condiction at outer boundary: phi(i=mpsi)=0
                endif
             endif
          else                      ! for dphi/dt solver
             if( i==0 ) this%gk_matid(it) = -1
             if( i==mpsi ) this%gk_matid(it) = -1
          endif
       end do
    end do
    if(id .ne. sum(mtheta)) call MPI_ABORT(mpi_comm_world,1,ierr)
    if(this%gk_matid(sum(mtheta) + mpsi + 2) .ne. -2) call MPI_ABORT(mpi_comm_world,1,ierr)
        
    ! form matrix data
    call poisson_initial( this )

  end subroutine new_phit_matrix
  
  !========================================================================
  !           subroutine delete_phit_matrix  
  !========================================================================
  subroutine delete_phit_matrix( this )
    implicit none
    type(phit_matrix_data) :: this

    deallocate(this%gk_matid)
    nullify(this%gk_matid)

  end subroutine delete_phit_matrix

  !========================================================================
  !          subroutine add_gk_value
  !========================================================================
  subroutine add_gk_value( this, ij0, ij, rs )
    use precision,only:wp,doubleprec
    use particle_decomp,only:myrank_partd,nproc_partd
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
    type(phit_matrix_data) :: this
    integer :: ij, ij0
    real(wp) rs 
    !
    integer nt, ierr, ii, jj
    PetscScalar val, zero

    ! insignificant point replaced by the original grid point !!!
    val = rs
    zero = 0.d0
    ii = this%gk_matid(ij0)
    jj = this%gk_matid(ij)
    if( ii >= 0 .and. jj >= 0 ) then ! BCs      
       if( mod(ii,nproc_partd) .eq. myrank_partd ) then
          if(val<0.00001)then
          !           if(val<1.0d-9)then
             call MatSetValues(this%data,1,ii,1,ii,val,ADD_VALUES,ierr)
             if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
             !ring(1,ij0,k)= this%ring(1,ij0,k)+rr
          else
             call MatSetValues(this%data,1,ii,1,jj,val,ADD_VALUES,ierr)
             if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
             ! symmmetrize for partitioning
             call MatSetValues(this%data,1,jj,1,ii,zero,ADD_VALUES,ierr)
             if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
          end if
       end if
    end if
 
  end subroutine add_gk_value

end module phit_matrix_class

!========================================================================
!            module solver_data_type
!========================================================================
!
module solver_data_type
  use precision,only:wp,doubleprec
  use phit_matrix_class
  implicit none
#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h" 
#include "finclude/petsckspdef.h"
#include "finclude/petscisdef.h"
  type solver_data
     Vec  xVec,bVec
     KSP  ksp
     VecScatter scat
     IS petsc_gts_is
     type(phit_matrix_data) :: mat
     Mat rhs_mat     
     !integer,dimension(:),pointer :: ipbc
  end type solver_data
end module solver_data_type
!
! solver_data_class
!
module solver_data_class
  use solver_data_type
  use particle_decomp,only:myrank_partd,nproc_partd,partd_comm
  implicit none
  public :: nullify_solver,destroy_solver,new_gk_solver
contains
  !========================================================================
  !               nullify_solver
  !========================================================================
  subroutine nullify_solver( this )
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h" 
#include "finclude/petscksp.h"
#include "finclude/petscis.h"
    type(solver_data) :: this    
    !
    this%mat%data = 0 !PETSC_NULL_OBJECT
    this%ksp = 0 !PETSC_NULL_OBJECT      
    this%scat = 0 !PETSC_NULL_OBJECT
    this%xVec = 0 !PETSC_NULL_OBJECT
    this%bVec = 0 !PETSC_NULL_OBJECT
    this%petsc_gts_is = 0 !PETSC_NULL_OBJECT
    this%rhs_mat = 0

    !nullify(this%ipbc)

  end subroutine nullify_solver
    !========================================================================
  !             subroutine solver_init
  !              - partition matrix    
  !              - create solver objects
  !========================================================================
  subroutine petsc_init()
    implicit none
!  #include "finclude/petsc.h"
#include "finclude/petscsys.h"
    integer ierror
    
    call PetscInitialize('./PETSc_input',ierror)
    
  end subroutine petsc_init
  !========================================================================
  !             subroutine solver_init
  !              - partition matrix    
  !              - create solver objects
  !========================================================================
  subroutine solver_init( this, prefix, Amat, neq, ierr )
    use precision,only: wp,doubleprec
    use global_parameters, only: mgrid, mype
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h" 
#include "finclude/petscksp.h"
#include "finclude/petscis.h"
#include "finclude/petscviewer.h"
    type(solver_data) :: this
    character,intent(in) :: prefix*(*)
    Mat :: Amat
    integer,intent(in) :: neq
    integer,intent(out) :: ierr
    IS               mis,isn,id_is,gts_petsc_is
    PetscInt         cnt(4096)
    PetscInt         npartitions
    Mat              Bmat
    KSP              ksp
    Vec              vec
    real(doubleprec):: tol,value(17),zeros(17)
    integer          i,flg,j,col(17),inz,n_loc,start,kk
    PetscViewer viewer
    MatPartitioning  mpart
    ! debug
    if( mype / nproc_partd .eq. -10 ) then
       call PetscViewerBinaryOpen(partd_comm, 'A0.m',FILE_MODE_WRITE,viewer,ierr)
       call MatView( Amat, viewer, ierr )
       call PetscViewerDestroy(viewer,ierr)
       write(*,*) " ---> debug matsum", ierr
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if( this%petsc_gts_is .eq. 0 ) then
       npartitions = nproc_partd
       call MatPartitioningCreate( partd_comm, mpart, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       call MatPartitioningSetAdjacency( mpart, Amat, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       call MatPartitioningSetFromOptions( mpart, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       call MatPartitioningApply( mpart, mis, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       call MatPartitioningDestroy( mpart, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       call ISPartitioningToNumbering( mis, gts_petsc_is, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
      ! call ISPartitioningCount( mis, cnt, ierr )
       call ISPartitioningCount( mis, npartitions, cnt, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       call ISDestroy(mis,ierr)
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )         
       ! side effect
       call ISInvertPermutation(gts_petsc_is,cnt(myrank_partd+1),this%petsc_gts_is,ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )         

       call ISSort( this%petsc_gts_is, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )

       call ISDestroy( gts_petsc_is, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
    else
       stop 'error 2'
    endif

    call ISAllGather( this%petsc_gts_is, isn, ierr )
    
!cj call MatGetSubMatrix(Amat,this%petsc_gts_is,isn,PETSC_DECIDE,MAT_INITIAL_MATRIX,Bmat,ierr)
    call MatGetSubMatrix(Amat,this%petsc_gts_is,this%petsc_gts_is,MAT_INITIAL_MATRIX,Bmat,ierr)
    if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )

    call MatDestroy( Amat, ierr )

    call ISDestroy( isn, ierr )

    Amat = Bmat

    if( this%bVec .eq. 0 ) then
       call MatGetOwnershipRange( Amat, start, kk, ierr )
       if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
       n_loc = kk - start
       
       call VecCreateMPI( partd_comm, n_loc, neq, this%xVec, ierr )
       call VecSetFromOptions( this%xVec, ierr )
       call VecDuplicate( this%xVec, this%bVec, ierr )
       
       call ISCreateStride(MPI_COMM_SELF,n_loc,start,1,id_is,ierr)
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       !     scatter object
       call VecCreateSeq( MPI_COMM_SELF, neq, vec, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       ! side effect 
       call VecScatterCreate(vec,this%petsc_gts_is,this%bVec,id_is,this%scat,ierr)
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       
       call ISDestroy( id_is, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
       call VecDestroy( vec, ierr )
       if(ierr .ne. 0) call MPI_ABORT( mpi_comm_world, 1, ierr )
    else
       stop 'SKIPPING vector allocation!!!'
    endif
    
    ! debug
    if( mype/nproc_partd .eq. -2 ) then
       call PetscViewerBinaryOpen(partd_comm, 'A.m',FILE_MODE_WRITE,viewer,ierr)
       call MatView( Amat, viewer, ierr )
       call PetscViewerDestroy(viewer,ierr)
       write(*,*) " ---> debug matsum", ierr
    endif
    ! KSP build
    call KSPCreate( partd_comm, ksp, ierr )
    call KSPSetOperators( ksp, Amat, Amat, DIFFERENT_NONZERO_PATTERN, ierr)
    tol = 1.d-5
    call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetType(ksp, KSPGMRES, ierr )
        
    call KSPSetOptionsPrefix( ksp, prefix, ierr )
    call KSPSetFromOptions( ksp, ierr )

    if( mype / nproc_partd .ne. 0 ) then
       call KSPMonitorCancel( ksp, ierr )
       !call KSPClearMonitor( ksp, ierr )
    endif    

    if(this%mat%data .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    this%mat%data = Amat
    this%ksp = ksp
    !
  end subroutine solver_init

  !========================================================================
  !               subroutine destroy_solver
  !========================================================================
  subroutine destroy_solver( this )
    implicit none
    type(solver_data) :: this
    integer ierr
    
    !  if( this%bVec .ne. 0 ) then
    call VecDestroy(this%xVec,ierr)
    
    !     if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecDestroy(this%bVec,ierr)
    !     if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    
    call VecScatterDestroy( this%scat, ierr )
    !     if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    
    call ISDestroy( this%petsc_gts_is, ierr )
    !  endif
    
    if( this%mat%data /= 0 ) then
       call MatDestroy( this%mat%data,ierr)
       !     if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    endif
    
    if( this%rhs_mat /= 0 ) then
       call MatDestroy( this%rhs_mat, ierr )
       !     if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    endif
    
    !  if( this%ksp .ne. 0 ) then
    call KSPDestroy(this%ksp,ierr)
    !     if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)            
    !  endif
    !
  end subroutine destroy_solver
  
  !========================================================================
  !              subroutine new_gk_solver
  !========================================================================
  subroutine new_gk_solver( this, prefix, mat_type )
    use precision,only:wp,doubleprec
    use global_parameters,only: mype,mtheta,mpsi,tite,ifelectron, &
                                lmin,lmax,rmin,rmax
    use field_array, only: psi_grid, igrid, gama0, gama1
    use phit_matrix_class
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscmat.h" 
#include "finclude/petscvec.h" 
    type(solver_data) :: this
    character*(*) prefix
    integer mat_type
    real(wp), external :: SMOOTHSTEP
    real(wp) instep,outstep
    ! locals
    integer iv,ib,sz,j,i
    type(phit_matrix_data) :: m
    Mat Amat
    integer ierr,neq,mindex,mring
    PetscScalar a
    ! init PETSc here, need to fix for multiple solvers!!!
    !call PetscInitialize('./.petscrc',ierr) 
    !if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    ! create matrix
    mring=2
    ! mindex=32*mring+1 ! estimate for matrix size
    mindex=64*mring+1 ! estimate for matrix size
    ! mindex=mtheta(mpsi)+64*mring+1 ! estimate for matrix size
    neq = sum( mtheta )    
  write(0,*)mype,' mring=',mring,' mindex=',mindex,' neq=',neq

    ! 
    !call MatCreateMPIAIJ(partd_comm,PETSC_DECIDE,PETSC_DECIDE,neq,neq,mindex,
    !   PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,Amat,ierr)
    call MatCreate( partd_comm, Amat, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)    
    call MatSetSizes(Amat, PETSC_DECIDE,PETSC_DECIDE, neq, neq, ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatSetFromOptions( Amat, ierr ) ! need ti separate for SuperLU
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    !call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex/5+1,PETSC_NULL_INTEGER,ierr)
    call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex/4+1,PETSC_NULL_INTEGER,ierr)
    call MatSeqAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    ! for GK Poisson equations (normalized):
    !     ETG case: (1+zeff*Te/Ti)*(Phi/Te)-tilde(Phi)/Te=delta-ne/n_e0
    !     ITG case: [1+(1/zeff)*Ti/Te]*(Phi/Ti)-tilde(Phi)/Ti=delta-ni/n_i0+(1/zeff)*<Phi>/Te
    if(mype .eq. 0) write(0,*) '[',mype,'] new_gk_solver: call integral form of PETSc operator, neq = ', neq
    call new_phit_matrix( m, Amat, mindex, mat_type )
    call delete_phit_matrix( m )
    
    ! I + tite - Phi_tilda
    iv = 0
    do i=0,mpsi
       sz = mtheta(i)
       write(0,*)mype,' mpsi:i=',i,' mtheta:i=',mtheta(i)

       !   ti/te     1
       if(mat_type==1)then
          if(ifelectron==0) then        ! adiabatic electron case
             a = -(tite(i) + 1.d0)
          elseif(ifelectron==1) then    ! full delta-f kinetic electron case
       !        a = -(0.01 + 1.d0)         ! full electron (no split)
       !        a = -(0.002 + 1.d0)         ! full electron (no split)
       !        a = -1.d0            

       !    ADDED Jan., 2010
              outstep=smoothstep(rmin,rmax,psi_grid(i))
       !       instep=smoothstep(lmin,lmax,psi_grid(i))
       !       instep=instep*(1.0d0-outstep)
       !       outstep=1.0d0-instep

       !       a = -(0.21 + 0.09*outstep + 1.d0)    ! full electron (no split)
                                                    ! --> a=-(0.21+1.0)  for r > rmax

              a = -(gama0 + gama1*outstep + 1.d0)    ! full electron (no split)







       !  ADDED TEST 
       !       a = -(0.01 + 1.d0)






       !       a = -(0.30 + 1.d0)                   ! full electron (no split)

       !    ADDED TEST Jan, 2010
       !       if(i.lt.144) then 
       !          a = -(0.01 + 1.d0)         ! full electron (no split)
       !       else
       !          a = -(0.10 + 1.d0)         ! added test Jan, 2010
       !       end if

          elseif(ifelectron==2) then    ! delta-h adiabatic electrons (split weight)
             a = -(tite(i) + 1.d0)   
          endif
       else         ! for dphi/dt solver

       !    a = -1.d0

       !    ADDED Jan., 2010
          outstep=smoothstep(rmin,rmax,psi_grid(i))
       !    a = -(0.21 + 0.09*outstep + 1.d0)    ! for solving dphi/dt in double-split weight
                                                 ! --> a=-(0.21+1.0)  for r > rmax

          a = -(gama0 + gama1*outstep + 1.d0)    ! for solving dphi/dt in double-split weight

       !  ADDED TEST 
       !       a = -(0.01 + 1.d0)
       endif
       do j=1,sz
          if(mod(i,nproc_partd) .eq. myrank_partd) then
             call MatSetValues(Amat,1,iv,1,iv,a,ADD_VALUES,ierr)
          endif
          iv = iv + 1
       enddo
       if(mype==0) then
          write(23,101) i, igrid(i), psi_grid(i), outstep, a
       end if
    enddo
101 format(2(2x, i8), 20(2x, ES13.6))
    if(iv .ne. neq) stop 'iv .ne. neq'

    ! communication phase
    call MatAssemblyBegin( Amat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyEnd( Amat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    !  A := 2 * I - A
    a = -1.d0
    call MatScale( Amat, a, ierr)
    call solver_init( this, prefix, Amat, neq, ierr )
  end subroutine new_gk_solver

!========================================================================
!         solve
!========================================================================
  subroutine solve( solver, nn, bb_gk, xx_gk, ierr )
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
    use global_parameters, only:mype,mtheta,mpsi,mgrid
    implicit none
    type(solver_data) :: solver
    integer,intent(in) :: nn
    integer,intent(out) :: ierr
    real(wp),intent(in) :: bb_gk(nn) 
    real(wp),intent(out) :: xx_gk(nn)
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq
    real(doubleprec),allocatable,dimension(:) :: xx, bb 
    integer count

    data count/0/
    count = count + 1

    neq = sum(mtheta)
    allocate(xx(neq),bb(neq))
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi 
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             bb(ivsl) = bb_gk(ivgts)
          elseif(bb_gk(ivgts).ne.bb(ivsl-itmax+1))then
          !            if(myrank_partd==0) print*, 'mype=', mype, ir, ivgts, ivsl, itmax, bb_gk(ivgts), bb(ivsl-itmax+1),ivsl-itmax+1, bb_gk(1),bb_gk(ivgts/2),bb(ivgts/2)
             stop 'bb_gk(ivgts).ne.bb(ivsl-itmax+1)'
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    
    ! solve GK Poisson equation
    call petsc_solve( solver, neq, bb, xx, ierr)

    ! set solution (output)
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts+1 
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             xx_gk(ivgts) = xx(ivsl)
          else
             xx_gk(ivgts) = xx(ivsl-itmax+1)
          endif
       enddo
    enddo

    deallocate(bb,xx) 

    ! debug
    if( count .eq. 1 .and. mype==-10 ) then
       call plot_array( mgrid, xx_gk, 'x.m' )
       call plot_array( mgrid, bb_gk, 'b.m' )
       stop 'plot'
    endif

  end subroutine solve

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
    integer :: nn,ierr
    real(doubleprec),intent(in) :: bb(nn)
    real(doubleprec),intent(out) :: xx(nn)
    type(solver_data) :: solver
    ! locals
    Vec              vec_yy
    integer          iplot,i,j
    PetscScalar t1
    PetscScalar,allocatable,dimension(:) :: yy, zz
    PetscViewer viewer

    if(solver%mat%data .eq. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    !     copy real vector into double precision
    allocate(yy(nn),zz(nn))
    yy = bb

    !     scatter rhs into PETSc vector
    call VecCreateSeqWithArray(MPI_COMM_SELF,1,nn,yy,vec_yy,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    
    !     call VecSet( solver%bVec, 1.e30, ierr ) ! debug
    call VecScatterBegin(solver%scat,vec_yy, solver%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(solver%scat,vec_yy, solver%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    
    !     Solve      
    call KSPSolve( solver%ksp, solver%bVec, solver%xVec, ierr )
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

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
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)      
    deallocate(yy,zz)
  end subroutine petsc_solve

end module solver_data_class

!========================================================================
!         subroutine plot_array( nn, arr )
!========================================================================
subroutine plot_array( nn, arr, str )
  use precision,only:wp,doubleprec
  implicit none
#include "finclude/petscsys.h"
#include "finclude/petscviewer.h"
  integer nn,ierr
  character :: str*(*)
  real(wp) arr(nn)
  PetscViewer viewer
  PetscScalar yy(nn)
  Vec v;
  
  yy = arr
  
  call VecCreateMPIWithArray(mpi_comm_self,nn,PETSC_DECIDE,yy,v,ierr)
  !
  call PetscViewerASCIIOpen(mpi_comm_self, str,viewer,ierr)
  call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
  call VecView( v, viewer, ierr )
  call PetscViewerDestroy(viewer,ierr)
  call VecDestroy( v, ierr )
  
end subroutine plot_array
