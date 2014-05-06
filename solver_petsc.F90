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

     Mat :: mdata ! mass matrix
     integer :: numbdy, bc_type
     integer,dimension(:),pointer :: bdy
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
                !   if( i==mpsi ) this%gk_matid(it) = -1         ! Dirichlet boundary condiction at outer boundary: phi(i=mpsi)=0
             else ! ZF_SEPARATION==0                
                if(ifelectron==0)  then ! adiabatic electrons
                   ! NOTE: any combination of following BCs produces a spike at the boundary for ad-electron case
                   ! if( i==0 ) this%gk_matid(it) = -1
                   ! if( i==mpsi ) this%gk_matid(it) = -1  ! Dirichlet boundary condiction at outer boundary: phi(i=mpsi)=0
                else
                   if( i==0 ) this%gk_matid(it) = -1
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
  use fem_data_class
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

     type(fem_data) :: fem

     !finite element
     ! order = 0 original solver
     ! order = 1 linear finite element
     ! order = 2 2nd order finite element
     ! order = 5 5th order c^1 finite element
     integer order
     real(doubleprec),allocatable,dimension(:) :: zero01, zero10
     Vec zeroVec01, zeroVec10
  end type solver_data
end module solver_data_type
!
! solver_data_class
!
module solver_data_class
  use solver_data_type
  use particle_decomp,only:myrank_partd,nproc_partd,partd_comm
  implicit none
  public :: nullify_solver,destroy_solver,new_gk_solver, new_fem_solver, new_fem_masssolver
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
    this%mat%mdata = 0 !PETSC_NULL_OBJECT
    this%ksp = 0 !PETSC_NULL_OBJECT      
    this%scat = 0 !PETSC_NULL_OBJECT
    this%xVec = 0 !PETSC_NULL_OBJECT
    this%bVec = 0 !PETSC_NULL_OBJECT
    this%petsc_gts_is = 0 !PETSC_NULL_OBJECT
    this%rhs_mat = 0
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
       write(*,*) " ---> debug matsum on mype", mype, nproc_partd
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
    call KSPSetOperators( ksp, Amat, Amat, SAME_PRECONDITIONER, ierr)
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
    
    if( this%mat%mdata /= 0 ) then
       call MatDestroy( this%mat%mdata,ierr)
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
    if(this%order .ge. 1) then
       call VecDestroy(this%zeroVec01,ierr)
       call VecDestroy(this%zeroVec10,ierr)
       deallocate(this%mat%bdy,this%mat%gk_matid)
    endif

  end subroutine destroy_solver
  
  !========================================================================
  !              subroutine new_gk_solver
  !========================================================================
  subroutine new_gk_solver( this, prefix, mat_type, bc_type )
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
    integer mat_type, bc_type
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
    !  write(0,*)mype,' mring=',mring,' mindex=',mindex,' neq=',neq

    ! 
    !call MatCreateMPIAIJ(partd_comm,PETSC_DECIDE,PETSC_DECIDE,neq,neq,mindex,
    !   PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,Amat,ierr)
    call MatCreate( partd_comm, Amat, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)    
    call MatSetSizes(Amat, PETSC_DECIDE,PETSC_DECIDE, neq, neq, ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatSetFromOptions( Amat, ierr ) ! need ti separate for SuperLU
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    !    call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex/5+1,PETSC_NULL_INTEGER,ierr)
    call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,ierr)
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
       !  write(0,*)mype,' mpsi:i=',i,' mtheta:i=',mtheta(i)

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
  !              subroutine new_fem_solver
  !========================================================================
  subroutine new_fem_solver( fem, this, prefix, mat_type, bc_type )
    use precision,only:wp,doubleprec
    use global_parameters, only: mype,mtheta,mpsi,numberpe
    use fem_data_class
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscmat.h" 
#include "finclude/petscvec.h" 
#include "finclude/petscviewer.h"
    type(fem_data) :: fem
    type(solver_data) :: this
    character,intent(in) :: prefix*(*)
    integer mat_type, bc_type
    ! locals
    real(wp) :: alpha,beta
    Mat Amat, Mmat
    integer ierr,neq,mindex,mring, i, j, ir, it, ib, ib2, id, ki
    PetscScalar a
    PetscViewer viewer
    integer, dimension(:), allocatable :: bdy
    Vec              vec_yy
    integer :: fe_ndof=1, ie, ind1, ind2, ind3, ind4, ind5, ind6, msize, nsize


    mring=2
    mindex=64*mring+1 ! estimate for matrix size
    neq = sum( mtheta )
    if(fem%order .gt. 2) fe_ndof=6
    neq = fem%nmesh_pts*fe_ndof
    if(mype==0) write(*,*) "	finite element: matrix dim neq=", sum(mtheta), neq

    call MatCreate( partd_comm, Amat, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)    
    call MatSetSizes(Amat, PETSC_DECIDE,PETSC_DECIDE, neq, neq, ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatSetFromOptions( Amat, ierr ) ! need ti separate for SuperLU
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    !    call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex/4+1,PETSC_NULL_INTEGER,ierr)
    call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,ierr)
    call MatSeqAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    call MatCreate( partd_comm, Mmat, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)    
    call MatSetSizes(Mmat, PETSC_DECIDE,PETSC_DECIDE, neq, neq, ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatSetFromOptions( Mmat, ierr ) ! need ti separate for SuperLU
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
!   call MatMPIAIJSetPreallocation(Mmat,mindex,PETSC_NULL_INTEGER,mindex/4+1,PETSC_NULL_INTEGER,ierr)
    call MatMPIAIJSetPreallocation(Mmat,mindex,PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,ierr)
    call MatSeqAIJSetPreallocation(Mmat,mindex,PETSC_NULL_INTEGER,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    alpha = 1.
    beta  = 1.
    !call form_tangent(Amat,Mmat,fem%iemax,fem%iet1,fem%iet2,&
    !     fem%iet3,neq,fem%points,alpha,beta,ierr)
    !call form_tangent(fem%order,Amat,Mmat,fem%iemax,fem%iet,&
    !     neq,fem%points,alpha,beta,ierr)
    call form_tangent(fem,Amat,Mmat,neq,alpha,beta,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    call MatAssemblyBegin( Amat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyEnd( Amat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyBegin( Mmat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyEnd( Mmat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    if(mype==0) write(*,*) "	finite elements matrix created bc_type=",bc_type

!cj reuse original code
    call new_gk_matid( this%mat, bc_type )

!cj imposing dirichlet bc
    ib=0
    if( bc_type .eq. 1 ) then
       ib = mtheta(0)
    else if( bc_type .eq. 2 ) then
       ib = mtheta(0) + mtheta(mpsi)
    else if( bc_type .eq. 3 ) then
       ib = mtheta(mpsi)
    endif
    if(fem%order .eq. 2) ib=ib*2
    this%mat%numbdy = ib*fe_ndof
    allocate(bdy(this%mat%numbdy))
    allocate(this%mat%bdy(this%mat%numbdy))
    it = 0
    ir = 0
    id = 0*fe_ndof
    ib2 = 0
    allocate(this%zero01(neq),this%zero10(neq))
    this%zero01=0.
    this%zero10=1.
    do id=1, sum( mtheta )
       if ( this%mat%gk_matid(id) .lt. 0 )  then
          ! use 0-based bdy index to set Dirichlet bc for lhs matrix
          do ki=0,fe_ndof-1
             ib2=ib2+1  !count bdy index and assign global index

             bdy(ib2)=(id-1)*fe_ndof+ki
             this%mat%bdy(ib2)=(id-1)*fe_ndof+ki

             ! use 1-based bdy index to set Dirichlet bc for rhs vector
             this%zero01((id-1)*fe_ndof+ki+1) = 1.
             this%zero10((id-1)*fe_ndof+ki+1) = 0.
          enddo
       endif
    end do
    if(mype==0) write(*,*) "	finite elements bdy ib ib2=", ib, ib2, mtheta(0), mtheta(mpsi)
    if(bc_type.ne.0 .and. fem%order .eq. 2) then
       !count the bdy vertex lying on triangle edges
       do ie=1,fem%iemax
          ind1= fem%iet(1,ie)
          ind2= fem%iet(3,ie)
          ind3= fem%iet(5,ie)

          if ( this%mat%gk_matid(ind1).lt.0 .and. this%mat%gk_matid(ind2).lt.0 ) then
             ib2=ib2+1
             ind4=fem%iet(2,ie)
             bdy(ib2)=ind4-1
             this%mat%bdy(ib2)=ind4-1
             this%zero01(ind4) = 1.
             this%zero10(ind4) = 0.
          endif
          if ( this%mat%gk_matid(ind2).lt.0 .and. this%mat%gk_matid(ind3).lt.0 ) then
             ib2=ib2+1
             ind5=fem%iet(4,ie)
             bdy(ib2)=ind5-1
             this%mat%bdy(ib2)=ind5-1
             this%zero01(ind5) = 1.
             this%zero10(ind5) = 0.
          endif
          if ( this%mat%gk_matid(ind3).lt.0 .and. this%mat%gk_matid(ind1).lt.0 ) then
             ib2=ib2+1
             ind6=fem%iet(6,ie)
             bdy(ib2)=ind6-1
             this%mat%bdy(ib2)=ind6-1
             this%zero01(ind6) = 1.
             this%zero10(ind6) = 0.
          endif
       enddo
       if(mype==0) write(*,*) "	finite elements o2 bdy ib ib2=", ib, ib2
    endif
    if(ib*fe_ndof .ne. ib2) call MPI_ABORT(mpi_comm_world,1,ierr)

!   Set Amat=1 on bdy vertex for applying Dirichlet bc
!   using zero-based bdy index
    if(ib .gt. 0) then
       a = 1.
       call MatZeroRows(Amat,ib2,bdy,a,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
       deallocate(bdy)
       if(mype==0) write(*,*) "	finite elements matrix bdy fix"
    endif

    call solver_init( this, prefix, Amat, neq,ierr)
    if(mype==0) write(0,*) '	finite elements solver created'
    
    ! finalize stiff matrix
    ! ... done in solver_init

    ! finalize mass matrix
    call MatGetSubMatrix(Mmat,this%petsc_gts_is,this%petsc_gts_is,MAT_INITIAL_MATRIX,this%mat%mdata,ierr)

    ! create zeroVec01 for applying bc
    call VecDuplicate(this%bVec, this%zeroVec01, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecCreateSeqWithArray(MPI_COMM_SELF,1,neq,this%zero01,vec_yy,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterBegin(this%scat,vec_yy, this%zeroVec01, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(this%scat,vec_yy, this%zeroVec01, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

    ! create zeroVec10 for applying bc
    call VecDuplicate(this%bVec, this%zeroVec10, ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecCreateSeqWithArray(MPI_COMM_SELF,1,neq,this%zero10,vec_yy,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterBegin(this%scat,vec_yy, this%zeroVec10, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(this%scat,vec_yy, this%zeroVec10, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    call VecDestroy( vec_yy, ierr )

    this%fem=fem
  end subroutine new_fem_solver


  !========================================================================
  !              subroutine new_fem_masssolver
  !========================================================================
  subroutine new_fem_masssolver( fem, this )
    use precision,only:wp,doubleprec
    use global_parameters, only: mype,mtheta,mpsi,numberpe
    use fem_data_class
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscmat.h" 
#include "finclude/petscvec.h" 
#include "finclude/petscviewer.h"
    type(fem_data) :: fem
    type(solver_data) :: this
    ! locals
    real(wp) :: alpha,beta
    Mat Amat, Mmat
    integer ierr,neq,mindex,mring


    mring=2
    mindex=64*mring+1 ! estimate for matrix size
    neq = sum( mtheta )
    if(fem%order .gt. 1) neq = fem%nmesh_pts
    if(mype==0) write(*,*) "	finite element: mass dim neq order=", sum(mtheta), neq, fem%order

    call MatCreate( partd_comm, Amat, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)    
    call MatSetSizes(Amat, PETSC_DECIDE,PETSC_DECIDE, neq, neq, ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatSetFromOptions( Amat, ierr ) ! need ti separate for SuperLU
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    !    call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex/4+1,PETSC_NULL_INTEGER,ierr)
    call MatMPIAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,ierr)
    call MatSeqAIJSetPreallocation(Amat,mindex,PETSC_NULL_INTEGER,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    call MatCreate( partd_comm, Mmat, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)    
    call MatSetSizes(Mmat, PETSC_DECIDE,PETSC_DECIDE, neq, neq, ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatSetFromOptions( Mmat, ierr ) ! need ti separate for SuperLU
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    !    call MatMPIAIJSetPreallocation(Mmat,mindex,PETSC_NULL_INTEGER,mindex/4+1,PETSC_NULL_INTEGER,ierr)
    call MatMPIAIJSetPreallocation(Mmat,mindex,PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,ierr)
    call MatSeqAIJSetPreallocation(Mmat,mindex,PETSC_NULL_INTEGER,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    alpha = 1.
    beta  = 1.
    !call form_tangent(Amat,Mmat,fem%iemax,fem%iet1,fem%iet2,&
    !     fem%iet3,neq,fem%points,alpha,beta,ierr)
    !call form_tangent(fem%order,Amat,Mmat,fem%iemax,fem%iet,&
    !     neq,fem%points,alpha,beta,ierr)
    call form_tangent(fem,Amat,Mmat,neq,alpha,beta,ierr)
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    call MatAssemblyBegin( Amat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyEnd( Amat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyBegin( Mmat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)
    call MatAssemblyEnd( Mmat, MAT_FINAL_ASSEMBLY, ierr )
    if(ierr .ne. 0) call MPI_ABORT(mpi_comm_world,1,ierr)

    call solver_init( this, ' ', Mmat, neq,ierr)
    if(mype==0) write(0,*) '	finite element: masssolver created'
    
    this%fem=fem

    !we don't need stiffness matrix
    call MatDestroy( Amat, ierr )
  end subroutine new_fem_masssolver


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
    use global_parameters, only:mype,mtheta,mpsi,mgrid,ZF_SEPARATION,ifelectron,ifdphi
    use field_array, only: phip00
    implicit none
    type(solver_data) :: solver
    integer,intent(in) :: nn
    integer,intent(out) :: ierr
    real(wp),intent(in) :: bb_gk(nn) 
    real(wp),intent(inout) :: xx_gk(nn)
    real(doubleprec),allocatable,dimension(:) :: xx_ho, bb_ho, hh, bh, hh_ho, bh_ho
    real(wp) :: neuma, neumt
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq
    real(doubleprec),allocatable,dimension(:) :: xx, bb 
    real(doubleprec),allocatable,dimension(:) :: bn,bn_ho
    integer :: count = 0


    count = count + 1
    neq = sum(mtheta)
    allocate(xx(neq),bb(neq))

    !set Neumann bc for rhs vector sept19-2013
    !at present it only works for order=1
    allocate(hh(neq),bh(neq),bn(neq))
    hh=0.
    bh=0.
    bn=0.

    ivsl = 0
    ivgts = 0
    do ir=0,mpsi 
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             bb(ivsl) = bb_gk(ivgts)
             xx(ivsl) = xx_gk(ivgts)
             hh(ivsl) = phip00(ir)
          elseif(bb_gk(ivgts).ne.bb(ivsl-itmax+1))then
             stop 'bb_gk(ivgts).ne.bb(ivsl-itmax+1)'
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'

    !cj2014.05.02 modified for solver_dphi_dt with homogeneous Dirichlet BC
    if(ifelectron==2.and.ifdphi==1) then
       xx(1:mtheta(0)) = 0.
       xx(neq-mtheta(mpsi)+1:neq) = 0.
    endif
    
    !cj2014.05.02 solver1 uses homogeneous Neumann BC, so don't call neumbcfix
#ifdef CJTEST
    if( solver%order .ge. 1) then
       if(solver%mat%bc_type==0 .or. solver%mat%bc_type==1) then
          it= sum( mtheta ) - mtheta(mpsi)
          ir=mpsi
          call neumbcfix(ir,mtheta(ir),it,neq,solver%fem,neuma,neumt,bn)
       endif
       if(solver%mat%bc_type==0 .or. solver%mat%bc_type==3) then
          it= 0
          ir=0
          call neumbcfix(ir,mtheta(ir),it,neq,solver%fem,neuma,neumt,bn)
       endif
    endif
#endif

    if( solver%fem%order .eq. 2) then
        ! interpolate linear to ho array
        allocate(xx_ho(solver%fem%nmesh_pts),bb_ho(solver%fem%nmesh_pts))
        call hogrid(solver%fem%iemax, solver%fem%iet, neq, xx, solver%fem%nmesh_pts, xx_ho)
        call hogrid(solver%fem%iemax, solver%fem%iet, neq, bb, solver%fem%nmesh_pts, bb_ho)

        ! cal bh_ho from hh(gts input) and add it to M*bb in subroutine petsc_solve
        allocate(hh_ho(solver%fem%nmesh_pts),bh_ho(solver%fem%nmesh_pts))
        call hogrid(solver%fem%iemax, solver%fem%iet, neq, hh, solver%fem%nmesh_pts, hh_ho)
        call add_hh2bb(solver%fem, solver%fem%nmesh_pts, hh_ho, bh_ho, ierr)
        ! set bh_ho to 0 if don't want to add it

        ! Neumann bc fix
        allocate(bn_ho(solver%fem%nmesh_pts))
        call hogrid(solver%fem%iemax, solver%fem%iet, neq, bn, solver%fem%nmesh_pts, bn_ho)
        bh_ho=bh_ho+bn_ho
        bn_ho=0.

    ! solve GK Poisson equation
    call petsc_solve( solver, solver%fem%nmesh_pts, bb_ho, xx_ho, bh_ho, ierr)

        ! interpolate linear from ho array
        call hogridinv(neq, xx, solver%fem%nmesh_pts, xx_ho)
        call hogridinv(neq, bb, solver%fem%nmesh_pts, bb_ho)
        deallocate(xx_ho,bb_ho)
        deallocate(hh_ho,bh_ho,bn_ho)

    else
        ! cal bh from hh(gts input) and add it to M*bb in subroutine petsc_solve
        if( solver%order.ge.1 .and. ZF_SEPARATION.eq.1) then
            call add_hh2bb(solver%fem, neq, hh, bh, ierr)
            !  set bh to 0 if don't want to add it
            !  bh=bh+bn
            !  bh=0.
        endif

    ! solve GK Poisson equation
    call petsc_solve( solver, neq, bb, xx, bh, ierr)
    endif

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

    deallocate(bb,xx,hh,bh,bn) 

    ! debug
    if( count .eq. 1 .and. mype==-10 ) then
       call plot_array( mgrid, xx_gk, 'x.m' )
       call plot_array( mgrid, bb_gk, 'b.m' )
       stop 'plot'
    endif

  end subroutine solve


!========================================================================
!         solve test phi-phitilda
!========================================================================
  subroutine phitilda_solve( solver, nn, bb_gk, xx_gk, ierr )
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
    use global_parameters, only:mype,mtheta,mpsi,mgrid,pi,zetamin,deltaz,deltat,R_0,a0,a1
    use field_array,only: psi_grid,igrid,qtinv
    use esi_interface, only: metric
    implicit none
    type(solver_data) :: solver
    integer,intent(in) :: nn
    integer,intent(out) :: ierr
    real(wp),intent(in) :: bb_gk(nn) 
    real(wp),intent(out) :: xx_gk(nn)
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq,ij
    real(doubleprec),allocatable,dimension(:) :: xx, bb
    real(doubleprec),allocatable,dimension(:) :: rhs, phi

    real(wp) zeta, pdum, dt, tdum, tmp
    real(wp), dimension(2) :: aaa,the,rj,zj,Bj
    real(wp) Mk, Lk, alfak
    real(wp) :: bfun1, bfun2, rhoj, Tj
    real(wp) zeta, pdum, dt, tdum, krho1, krho2, fip, tmp1, tmp2
    real(wp), external :: temperature
    integer, external :: metric_elements
    real(wp),external :: besei0
    integer metr

    integer :: count = 0


    count = count + 1
    neq = sum(mtheta)
    allocate(xx(neq),bb(neq))
    allocate(rhs(nn), phi(nn))
    xx=0.0_wp
    ivsl = 0
    ivgts = 0
    zeta = zetamin + deltaz ! zeta for my plane   
    do ir=0,mpsi 
       itmax=mtheta(ir)
       pdum=psi_grid(ir)
       dt = deltat(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          ij = igrid(ir) + it

             tdum = real(it,wp)*dt + zeta*qtinv(ir)
             tdum=tdum/(2.0_wp*pi)
             tdum=2.0_wp*pi*(tdum-aint(tdum))

             Mk=10.
             Lk=10.
             alfak=0.
             phi(ij) = sin(pi*Mk*(pdum-a0)/(a1-a0))*cos(Lk*tdum+alfak)

             aaa(1)=pdum
             the(1)=tdum
             call esigetrzb(rj,zj,Bj,aaa,the,1)

             Tj = temperature( pdum, 0.0_wp )
             metr=metric_elements(pdum,tdum,zeta)
             rhoj=Tj/Bj(1)/Bj(1)
             fip=pi*Mk/(a1-a0)

             ! exp(-tmp) multiplied bessel function (k_perp^2*rhoj^2) for mode 1
             tmp1 = fip*fip * metric(1) + Lk * Lk * metric(2) + 2.0d0 * fip * Lk * metric(4)
             krho1 = tmp1 * rhoj
             bfun1 = besei0( krho1 )

             ! e myltiplied bessel function (k_perp^2*rhoj^2) for mode 2
             tmp2 = fip *fip * metric(1) + Lk * Lk * metric(2) - 2.d0 * fip * Lk * metric(4)
             krho2 = tmp2 * rhoj
             bfun2 = besei0( krho2 )

             ! phi - phi_tilda
             rhs(ij) = sin( fip*(pdum-a0) + Lk*tdum+alfak ) * (1. - bfun1) /2. + &
                       sin( fip*(pdum-a0) - Lk*tdum-alfak ) * (1. - bfun2) /2.

             if(ir==0 .or. ir==mpsi) rhs(ij) = phi(ij)
       enddo
           rhs(igrid(ir))=rhs(igrid(ir)+mtheta(ir)) !(:,igrid(i)) is missing in last do-loop
    enddo

    neq = sum(mtheta)
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi 
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             bb(ivsl) = rhs(ivgts)
          elseif(rhs(ivgts).ne.bb(ivsl-itmax+1))then
             stop 'bb_gk(ivgts).ne.bb(ivsl-itmax+1)'
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    
    ! solve GK Poisson equation
    call petsc_solve( solver, neq, bb, xx, bb, ierr)

    ! set solution (output)
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi
       itmax=mtheta(ir)
       pdum=psi_grid(ir)
       dt = deltat(ir)
       do it=0,itmax
          ivgts = ivgts+1 
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             xx_gk(ivgts) = xx(ivsl)

             tdum = real(it,wp)*dt + zeta*qtinv(ir)
             tdum=tdum/(2.0_wp*pi)
             tdum=2.0_wp*pi*(tdum-aint(tdum))
             tmp = sin(pi*Mk*(pdum-a0)/(a1-a0))*cos(Lk*tdum+alfak)

             if( mype.eq.0 .and. (ir.ne.0 .and. ir.ne.mpsi)) &
             write(111,111) ir, it, ivsl, xx(ivsl), tmp
111 format("finite element test solve=", 3(2x,i7), 2(1x,d12.6))

          else
             xx_gk(ivgts) = xx(ivsl-itmax+1)
          endif
       enddo
    enddo

    deallocate(rhs,xx,bb,phi) 
    call PetscFinalize( ierr )
    stop

    ! debug
    if( count .eq. 1 .and. mype==-10 ) then
       call plot_array( mgrid, xx_gk, 'x.m' )
       call plot_array( mgrid, bb_gk, 'b.m' )
       stop 'plot'
    endif
  end subroutine phitilda_solve


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
    use global_parameters, only:mype,mtheta,mpsi,mgrid,zetamin,deltaz,deltat,R_0,a0,a1,pi,csc0,m0ms,c0cs_i,msm0_i
    use field_array,only: psi_grid,qtinv
    implicit none
#include "finclude/petscsys.h"
    type(solver_data) :: solver
    integer :: nn
    integer :: ierr
    ! locals
    integer :: ir,it,ivgts,ivsl,ivsl2,itq,itmax,ie,iv,jv,kv,neq,mt
    real(doubleprec),allocatable,dimension(:) :: xx, bb, hh, bh, xx1
    real(doubleprec),allocatable,dimension(:) :: xx_ho, bb_ho, hh_ho, bh_ho
    real(doubleprec),allocatable,dimension(:) :: bn,bn_ho
    real(wp) :: neuma, neumt
    real(wp) zeta, pdum, dt, tdum, tmp
    real(wp), dimension(2) :: aaa,the,rj,zj,Bj
    real(wp) cf1, cf2, cf3, cf4, cf5
    real(wp) alpha, beta, ns0, ne0, gs, vthe3, Zs
    real(wp), external ::temperature_e, density
    integer, external :: metric_elements

    real(wp) M, L, alfa
    real(wp) :: rms=0., eps, err, pim, rr, zz
    logical :: flg
    integer :: whichone = 5
    integer :: ki, fe_ndof, i0, i1, i2, i3, i4, i5


    if(mype.eq.0) &
    write(0, *) "	finite element test solve", whichone, "of five"
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-whichone',whichone,flg,ierr)

    fe_ndof=1
    if( solver%fem%order .gt. 2) fe_ndof=6

!cj the following is a combianation of solve and tri_mesh
    neq = sum(mtheta)
    allocate(xx(neq*fe_ndof),bb(neq*fe_ndof))
    allocate(xx1(neq*fe_ndof))
    allocate(hh(neq*fe_ndof), bh(neq*fe_ndof), bn(neq*fe_ndof))
    xx=0.
    hh=1.
    bh=0.
    bn=0.  !dphidt related
    xx1=0.
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

             M=10.
             L=10.
             alfa=0.

             tdum = real(it)*dt + zeta*qtinv(ir)
             tdum=tdum/(2.0_wp*pi)
             tdum=2.0_wp*pi*(tdum-aint(tdum))

             aaa(1)=pdum
             the(1)=tdum
             call esigetrzb(rj,zj,Bj,aaa,the,1)

                i0=(ivsl-1)*fe_ndof+1
                if(fe_ndof .gt. 1) then
                   i1=(ivsl-1)*fe_ndof+2
                   i2=(ivsl-1)*fe_ndof+3
                   i3=(ivsl-1)*fe_ndof+4
                   i4=(ivsl-1)*fe_ndof+5
                   i5=(ivsl-1)*fe_ndof+6
                endif

             if(whichone==1) then
                !u = (x + y)^2,
                !rhs: c1f*d^2/dx^2 + cf2*d^/dxdy + cf3*d^2/dy^2 + cf4*d/dx + cf5*d/dy = 4
                cf1=1.
                cf2=0.
                cf3=1.
                xx(i0) = ( (rj(1)-R_0) + zj(1) ) * ( (rj(1)-R_0) + zj(1) )
                bb(i0) = cf1*2. + cf2*2.*2. + cf3*2.
                if(fe_ndof .gt. 1) then
                   xx(i1) = ( (rj(1)-R_0) + zj(1) ) * 2.
                   xx(i2) = ( (rj(1)-R_0) + zj(1) ) * 2.
                   xx(i3) = 2.
                   xx(i4) = 2.
                   xx(i5) = 2.
                   bb(i1) = 0.
                   bb(i2) = 0.
                   bb(i3) = 0.
                   bb(i4) = 0.
                   bb(i5) = 0.
                endif
             else if(whichone==2) then
                !u = sin(pi M x) cos( L y), del u = -( (pi M)^2 + L^2 ) u
                pim = pi*M
                rr=rj(1)-R_0
                zz=zj(1)
                xx(i0) = sin(pim*rr) * cos(L*zz+alfa)
                bb(i0) = (-pim*pim - L*L) * xx(i0) 
                if(fe_ndof .gt. 1) then
                   xx(i1) = cos(pim*rr) * cos(L*zz+alfa) * (pim)
                   xx(i2) = sin(pim*rr) * sin(L*zz+alfa) * (-L)
                   xx(i3) = -(pim*pim) * xx(i0)
                   xx(i4) = cos(pim*rr) * sin(L*zz+alfa) * (-pim*L)
                   xx(i5) = -(L*L) * xx(i0)
                   bb(i1) = (-(pim)*(pim) - L*L) * xx(i1) 
                   bb(i2) = (-(pim)*(pim) - L*L) * xx(i2) 
                   bb(i3) = (-(pim)*(pim) - L*L) * xx(i3) 
                   bb(i4) = (-(pim)*(pim) - L*L) * xx(i4) 
                   bb(i5) = (-(pim)*(pim) - L*L) * xx(i5) 
                endif
             else if(whichone==3) then
                !in at coordinates
                !matrix: Amatrix = dnx[i]*dnx[j] + dny[i]*dny[j]
                !      d^2 phi   d^2 phi    d^2 phi
                !del u=------- + ------- + --------
                !        da^2      dt^2      da dt

                !u = sin(pi M (a-a0)/(a1-a0)) cos( L theta + alfa),
                !rhs: del u = -( (pi M/(a1-a0))^2 + L^2 )u - (pi M/(a1-a0) L)cos sin
                pim = pi*M/(a1-a0)
                rr=pdum-a0
                zz=tdum
                xx(i0) = sin(pim*rr) * cos(L*zz+alfa)
                bb(i0) = (-pim*pim - L*L) * xx(i0) - (2.*pim*L)*cos(pim*rr)*sin(L*zz+alfa)
                if(fe_ndof .gt. 1) then
                   xx(i1) = cos(pim*rr) * cos(L*zz+alfa) * (pim)
                   xx(i2) = sin(pim*rr) * sin(L*zz+alfa) * (-L)
                   xx(i3) = -(pim*pim) * xx(i0)
                   xx(i4) = cos(pim*rr) * sin(L*zz+alfa) * (-pim*L)
                   xx(i5) = -(L*L) * xx(i0)
                   bb(i1) = (-(pim)*(pim) - L*L) * xx(i1) - (2.*pim*L)*sin(pim*rr)*sin(L*zz+alfa) * (-pim)
                   bb(i2) = (-(pim)*(pim) - L*L) * xx(i2) - (2.*pim*L)*cos(pim*rr)*cos(L*zz+alfa) * (L)
                   bb(i3) = (-(pim)*(pim) - L*L) * xx(i3) - (2.*pim*L)*cos(pim*rr)*sin(L*zz+alfa) * (-pim*pim)
                   bb(i4) = (-(pim)*(pim) - L*L) * xx(i4) - (2.*pim*L)*sin(pim*rr)*cos(L*zz+alfa) * (-pim*L)
                   bb(i5) = (-(pim)*(pim) - L*L) * xx(i5) - (2.*pim*L)*cos(pim*rr)*sin(L*zz+alfa) * (-L*L)
                endif

             else if(whichone==4) then
                !in at coordinates
                !matrix: Amatrix = 
                !      d^2 phi / da ^2  da ^2\     d^2 phi / da  dt     da  dt \   d^2 phi / dt ^2  dt ^2\
                !del u=------- |(--)  +(--)  | + 2 ------- |(--)(--) + (--)(--)| + ------- |(--)  +(--)  |
                !      da^2    \ dx     dy   /     da dt   \ dx  dx     dy  dy /   dt^2    \ dx     dy   /

                xx(i0) = pdum*pdum
                bb(i0) = 4.
                if(fe_ndof .gt. 1) then
                   xx(i1) = pdum*2.
                   xx(i2) = 0.
                   xx(i3) = 2.
                   xx(i4) = 0.
                   xx(i5) = 0.
                   bb(i1) = 0.
                   bb(i2) = 0.
                   bb(i3) = 0.
                   bb(i4) = 0.
                   bb(i5) = 0.
                endif
             else
                !u = sin(pi M (a-a0)/(a1-a0)) cos( L theta + alfa),
                !rhs: del u = -( (pi M / (a1-a0))^2 + L^2 ) u
                xx(ivsl) = sin(pi*M*(pdum-a0)/(a1-a0))*cos(L*tdum+alfa)

                vthe3=temperature_e(pdum,0.0_wp)
                alpha = abs(csc0) * m0ms / vthe3
                bb(ivsl) = alpha * xx(ivsl)

                ns0=density(pdum,0.0_wp)
                ne0=density(pdum,0.0_wp)
                beta = -1. / ne0
                Zs=1.
                gs = Zs / Bj(1) / Bj(1) * ns0 * c0cs_i * msm0_i
                bb(ivsl) = 1.
             endif

             ! assign the boundary vertex
             ! must be set to negative xx since Amat was set to -1 on bdy vertex
             ! if( ir==0 .or. ir==mpsi ) then
             do ki=1,fe_ndof
                ivsl2=(ivsl-1)*fe_ndof+ki
                if( solver%zero01(ivsl2) .gt. 0 ) then
                   xx(ivsl2) =   xx(ivsl2)
                   !if(mype.eq.0) write(*,119) ir,it,ivsl,ivsl2, xx(ivsl2)
                   !119 format("o6=", 4(2x,i7), (3x,d10.4))
                endif
             enddo
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'

    !set Neumann bc for rhs vector sept19-2013
    !at present it only works for order=1
    if(solver%mat%bc_type==0 .or. solver%mat%bc_type==1) then
       it= sum( mtheta ) - mtheta(mpsi)
       ir=mpsi
       call neumbcfix(ir,mtheta(ir),it,neq,solver%fem,neuma,neumt,bn)
    endif
    if(solver%mat%bc_type==0 .or. solver%mat%bc_type==3) then
       it= 0
       ir=0
       call neumbcfix(ir,mtheta(ir),it,neq,solver%fem,neuma,neumt,bn)
    endif

    if( solver%fem%order .eq. 2) then
        if(mype==0) &
        write(*,*) '	order neq nmesh_pts=', &
                   solver%fem%order, neq, solver%fem%nmesh_pts

        allocate(xx_ho(solver%fem%nmesh_pts),bb_ho(solver%fem%nmesh_pts))
        allocate(hh_ho(solver%fem%nmesh_pts),bh_ho(solver%fem%nmesh_pts))
        call hogrid(solver%fem%iemax, solver%fem%iet, neq, xx, solver%fem%nmesh_pts, xx_ho)
        call hogrid(solver%fem%iemax, solver%fem%iet, neq, bb, solver%fem%nmesh_pts, bb_ho)

        ! cal bh_ho from hh(gts input) and add it to M*bb in subroutine petsc_solve
        call hogrid(solver%fem%iemax, solver%fem%iet, neq, hh, solver%fem%nmesh_pts, hh_ho)
        call add_hh2bb(solver%fem, solver%fem%nmesh_pts, hh_ho, bh_ho, ierr)
        ! set bh_ho to 0 if don't want to add it

        ! Neumann bc fix
        allocate(bn_ho(solver%fem%nmesh_pts))
        call hogrid(solver%fem%iemax, solver%fem%iet, neq, bn, solver%fem%nmesh_pts, bn_ho)
        bh_ho=bh_ho+bn_ho
        bh_ho=0.

    ! solve GK Poisson equation
    call petsc_solve( solver, solver%fem%nmesh_pts, bb_ho, xx_ho, bh_ho, ierr)

        call hogridinv(neq, xx, solver%fem%nmesh_pts, xx_ho)
        call hogridinv(neq, bb, solver%fem%nmesh_pts, bb_ho)

        deallocate(xx_ho,bb_ho,hh_ho,bh_ho)

    else
        ! cal bh from hh(gts input) and add it to M*bb in subroutine petsc_solve
        call add_hh2bb(solver%fem, neq*fe_ndof, hh, bh, ierr)
        ! set bh to 0 if don't want to add it

        ! Neumann bc fix
        bh=bh+bn
        bh=0.

    call petsc_solve( solver, neq*fe_ndof, bb, xx, bh, ierr)
    endif

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
             tdum=tdum/(2.0_wp*pi)
             tdum=2.0_wp*pi*(tdum-aint(tdum))

             aaa(1)=pdum
             the(1)=tdum
             call esigetrzb(rj,zj,Bj,aaa,the,1)

                i0=(ivsl-1)*fe_ndof+1
                if(fe_ndof .gt. 1) then
                   i1=(ivsl-1)*fe_ndof+2
                   i2=(ivsl-1)*fe_ndof+3
                   i3=(ivsl-1)*fe_ndof+4
                   i4=(ivsl-1)*fe_ndof+5
                   i5=(ivsl-1)*fe_ndof+6
                endif
             if(whichone==1) then
                xx1(i0) = ( (rj(1)-R_0) + zj(1) ) * ( (rj(1)-R_0) + zj(1) )
                if(fe_ndof .gt. 1) then
                   xx1(i1) = ( (rj(1)-R_0) + zj(1) ) * 2.
                   xx1(i2) = ( (rj(1)-R_0) + zj(1) ) * 2.
                   xx1(i3) = 2.
                   xx1(i4) = 2.
                   xx1(i5) = 2.
                endif
             else if(whichone==2) then
                xx1(i0) = sin(pi*M*(rj(1)-R_0)) * cos(L*zj(1)+alfa)
                if(fe_ndof .gt. 1) then
                   xx1(i1) = cos(pi*M*(rj(1)-R_0)) * cos(L*zj(1)+alfa) * (pi*M)
                   xx1(i2) = sin(pi*M*(rj(1)-R_0)) * sin(L*zj(1)+alfa) * (-L)
                   xx1(i3) =-(pi*M)*(pi*M) * xx1(i0)
                   xx1(i4) = cos(pi*M*(rj(1)-R_0)) * sin(L*zj(1)+alfa) * (-pi*M*L)
                   xx1(i5) = -(L*L) * xx1(i0)
                endif
             else if(whichone==3) then
                pim = pi*M/(a1-a0)
                xx1(i0) = sin(pim*(pdum-a0))*cos(L*tdum+alfa)
                if(fe_ndof .gt. 1) then
                   xx1(i1) = cos(pim*(pdum-a0)) * cos(L*tdum+alfa) * (pim)
                   xx1(i2) = sin(pim*(pdum-a0)) * sin(L*tdum+alfa) * (-L)
                   xx1(i3) = -(pim*pim) * xx(i0)
                   xx1(i4) = cos(pim*(pdum-a0)) * sin(L*tdum+alfa) * (-pim*L)
                   xx1(i5) = -(L*L) * xx(i0)
                endif
             else if(whichone==4) then
                xx1(i0) = pdum*pdum
                if(fe_ndof .gt. 1) then
                   xx1(i1) = pdum*2.
                   xx1(i2) = 0.
                   xx1(i3) = 2.
                   xx1(i4) = 0.
                   xx1(i5) = 0.
                endif
             else
                tmp = sin(pi*M*(pdum-a0)/(a1-a0))*cos(L*tdum+alfa)
             endif

             ! print the inner vertex for debugging
             rms = rms &
                   +abs(xx1(i0)-xx(i0))/abs(xx1(i0)+xx(i0))*abs(xx1(i0)-xx(i0))/abs(xx1(i0)+xx(i0))
                if(fe_ndof .gt. 1) then
             rms = rms &
                   +abs(xx1(i1)-xx(i1))/abs(xx1(i1)+xx(i1))*abs(xx1(i1)-xx(i1))/abs(xx1(i1)+xx(i1)) &
                   +abs(xx1(i2)-xx(i2))/abs(xx1(i2)+xx(i2))*abs(xx1(i2)-xx(i2))/abs(xx1(i2)+xx(i2)) &
                   +abs(xx1(i3)-xx(i3))/abs(xx1(i3)+xx(i3))*abs(xx1(i3)-xx(i3))/abs(xx1(i3)+xx(i3)) &
                   +abs(xx1(i4)-xx(i4))/abs(xx1(i4)+xx(i4))*abs(xx1(i4)-xx(i4))/abs(xx1(i4)+xx(i4)) &
                   +abs(xx1(i5)-xx(i5))/abs(xx1(i5)+xx(i5))*abs(xx1(i5)-xx(i5))/abs(xx1(i5)+xx(i5))
                endif
             err= abs(xx1(i0)-xx(i0))/abs(xx1(i0)+xx(i0))
                if(fe_ndof .gt. 1) then
             err= err &
                 +abs(xx1(i1)-xx(i1))/abs(xx1(i1)+xx(i1)) &
                 +abs(xx1(i2)-xx(i2))/abs(xx1(i2)+xx(i2)) &
                 +abs(xx1(i3)-xx(i3))/abs(xx1(i3)+xx(i3)) &
                 +abs(xx1(i4)-xx(i4))/abs(xx1(i4)+xx(i4)) &
                 +abs(xx1(i5)-xx(i5))/abs(xx1(i5)+xx(i5))
                endif
             eps=.1
             if( mype.eq.0 ) then
!                if(fe_ndof .gt. 1) then
!             write(111,111) ir, ivsl, xx1(i0), xx(i0),  &
!                                      xx1(i1), xx(i1),  &
!                                      xx1(i2), xx(i2),  &
!                                      xx1(i3), xx(i3),  &
!                                      xx1(i4), xx(i4),  &
!                                      xx1(i5), xx(i5), err!/fe_ndof
!111 format("c1=", 2(1x,i6), 13(2x,d10.4))
!                else
!             write(111,112) ir, ivsl, xx1(i0), xx(i0), err!/fe_ndof
!112 format("c1=", 2(1x,i6),  3(2x,d10.4))
!                endif
                write(111,111) ivsl, xx1(i0), xx(ivsl), pdum, tdum, ir, it, itmax, err
111 format("finite element test solve=", 1(2x,i7), 4(2x,d10.4), 3(2x,i5), (3x,d10.4))
             endif
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    call flush(111)
    if(mype.eq.0) write(0, *) "	finite element test rms=",sqrt(rms)/ivsl
    if(mype.eq.0) write(0, *) "	finite element test solve end"
    
    deallocate(bb,xx) 
    deallocate(xx1,hh,bh,bn) 

    ! debug
    if( whichone .eq. 0 .and. mype==-10 ) then
       call plot_array( neq, xx, 'x.m' )
       call plot_array( neq, bb, 'b.m' )
       stop 'plot'
    endif

  end subroutine testsolve


!========================================================================
!         testmasssolve
!========================================================================
  subroutine testmasssolve( masssolver )
    !     solve M (dxdr) = u y Pijk
    !           M (dxdz) =-u x Pijk

    !use global_parameters
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mtheta,mpsi,mgrid,pi,zetamin,deltaz,deltat,R_0,a0,a1
    use field_array,only: psi_grid,qtinv
    implicit none
    type(solver_data) :: masssolver
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,mt
    real(doubleprec),allocatable,dimension(:) :: xx, bb, fcj, dxdr, dxdz
    real(doubleprec),allocatable,dimension(:) :: xx_ho, bb_ho, fcj_ho
    real(wp) zeta, pdum, dt, tdum, tmp
    real(wp),dimension(2) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
    integer, external :: metric_elements
    integer :: ierr

    real(wp) :: M, L, alfa, rms=0., pim, rr, zz, tmpr, tmpz
    integer :: i,ii,j,jj,k,kk
    integer :: fe_ns=3


    if(masssolver%order .gt. 0) then
       if(mype.eq.0) write(0, *) "	mass solver:: not implemented for high order elements"
       return
    else
       if(mype.eq.0) write(0, *) "	mass solver:: find dxdr dxdz begin"
    endif

!cj the following is a combianation of solve and tri_mesh
    allocate(xx(mgrid),bb(mgrid))
    allocate(dxdr(mgrid),dxdz(mgrid),fcj(mgrid))
    bb=0.
    dxdr=0.
    dxdz=0.
    ivsl = 0
    ivgts = 0
    zeta = zetamin + deltaz ! zeta for my plane   
    do ir=0,mpsi 
       itmax=mtheta(ir)
       pdum=psi_grid(ir)
       dt = deltat(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          tdum = real(it)*dt + zeta*qtinv(ir)
          tdum=tdum/(2.0_wp*pi)
          tdum=2.0_wp*pi*(tdum-aint(tdum))

          aaa(1)=pdum
          the(1)=tdum
          call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)

         !whichone = 1
         !bb(ivgts) = 1.
          bb(ivgts) =   (rj(1)-R_0) + zj(1)  
         !bb(ivgts) = ( (rj(1)-R_0) + zj(1) ) * ( (rj(1)-R_0) + zj(1) )
         !fcj(ivgts) = 1.

         !whichone = 3
         !bb(ivgts) = pdum*pdum
          M=1
          L=1
          alfa=1
          pim = pi*M/(a1-a0)
          rr=pdum-a0
          zz=tdum
         !bb(ivgts) = sin(pim*rr) * cos(L*zz+alfa)
         !fcj(ivgts) = (zaj(1)*rqj(1)-zqj(1)*raj(1))/rj(1)
       enddo
    enddo
    if(ivgts.ne.mgrid) then
       if(mype.eq.0) write(*,*) "mgrid=",mgrid,"ivgts=",ivgts
       stop'mass solver: ivgts.ne.mgrid'
    endif

    call find_dxdadxdt( masssolver, mgrid, bb, dxdr, dxdz, ierr )
    !call find_dxda( masssolver, mgrid, bb, dxdr, ierr )
    !call find_dxdt( masssolver, mgrid, bb, dxdz, ierr )

    ivsl = 0
    ivgts = 0
    zeta = zetamin + deltaz ! zeta for my plane   
    do ir=0,mpsi 
       itmax=mtheta(ir)
       pdum=psi_grid(ir)
       dt = deltat(ir)
       do it=0,itmax
          ivgts = ivgts + 1

          tdum = real(it)*dt + zeta*qtinv(ir)
          tdum=tdum/(2.0_wp*pi)
          tdum=2.0_wp*pi*(tdum-aint(tdum))

          aaa(1)=pdum
          the(1)=tdum
          call esigetrzb(rj,zj,Bj,aaa,the,1)

         !tmp = 0.
          tmpr = 1.
          tmpz = 1.
         !tmpr = 2. * ( (rj(1)-R_0) + zj(1) )
         !tmpz = 2. * ( (rj(1)-R_0) + zj(1) )
         !tmp = 2.*pdum
          M=1
          L=1
          alfa=1
          pim = pi*M/(a1-a0)
          rr=pdum-a0
          zz=tdum
         !tmpr = cos(pim*rr) * cos(L*zz+alfa) * (pim)
         !tmpz = sin(pim*rr) * sin(L*zz+alfa) * (-L)

          if(mype.eq.0)write(115,115) ivgts, bb(ivgts), tmpr, dxdr(ivgts), tmpz, dxdz(ivgts)
          115 format("mass solver=", 1(2x,i7), 5(1x,d12.6))
       enddo
    enddo
    if(ivgts.ne.mgrid)stop'mass solver: ivgts.ne.mgrid'
    call flush(115)
    if(mype.eq.0) write(0, *) &
      "	mass solver: rms(dxdr dxdz)=",sqrt(sum(dxdr*dxdr))/mgrid, sqrt(sum(dxdz*dxdz))/mgrid
    
    deallocate(xx,bb) 
    deallocate(dxdr,dxdz,fcj)

  end subroutine testmasssolve


!========================================================================
!         testpetscoption
!========================================================================
  subroutine testAxb( solver )
    !     solve Ax=b with x=1, b=Ax
    !     
    !use global_parameters
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mtheta
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
    type(solver_data) :: solver
    integer :: nn, ierr
    ! locals
    real(doubleprec),allocatable,dimension(:) :: xx, bb
    real :: neq
    logical :: flg
    real(wp) norm
    Vec xVec, bVec;
    PetscInt its; 
    KSPConvergedReason reason;


    if(mype.eq.0) write(0, *) "	testAxb:"
    neq = sum(mtheta)

    call VecSet( solver%xVec, 1.0, ierr )
    call VecDuplicate( solver%xVec, xVec, ierr )
    call VecCopy( solver%xVec, xVec, ierr )
        call VecNorm(xVec,NORM_1,norm, ierr);
        if(mype.eq.0) write(0, *) "	testAxb x=",norm/neq

    call VecSet( solver%bVec, 0.0, ierr )
    call MatMult(solver%mat%data, solver%xVec, solver%bVec, ierr)
        call VecNorm(solver%bVec,NORM_1,norm, ierr);
        if(mype.eq.0) write(0, *) "	testAxb b=",norm/neq

    call VecSet( solver%xVec, 0.0, ierr )
    call KSPSolve( solver%ksp, solver%bVec, solver%xVec, ierr )

        call KSPGetIterationNumber(solver%ksp,its,ierr);CHKERRQ(ierr);
        if(mype.eq.0) write(0, *) "	testAxb its=",its

        call VecNorm(solver%xVec,NORM_1,norm, ierr);
        if(mype.eq.0) write(0, *) "	testAxb x=",norm/neq

        call KSPGetConvergedReason(solver%ksp,reason,ierr);
        if(mype.eq.0) write(0, *) "	testAxb reason=",reason

        norm=0.
        call VecDuplicate( solver%bVec, bVec, ierr )
        call MatMult(solver%mat%data,solver%xVec,bVec,ierr);
        call VecNorm(bVec,NORM_1,norm, ierr);
        if(mype.eq.0) write(0, *) "	testAxb b=",norm/neq
        call VecAXPY(bVec,-1.0,solver%bVec,ierr);
        call VecNorm(bVec,NORM_2,norm,ierr);
        if(mype.eq.0) write(0, *) "	testAxb residual norm=",norm/sqrt(neq)

    call VecDestroy( xVec, ierr )
    call VecDestroy( bVec, ierr )
    call VecSet(solver%xVec, 0., ierr)
    call VecSet(solver%bVec, 0., ierr)

    if(mype.eq.0) write(0, *) "	testAxb end"
    
  end subroutine testAxb


!========================================================================
!         find dxdr dxdz
!========================================================================
  subroutine find_dxdadxdt( solver, nn, bb_gk, xa_gk, xt_gk, ierr )

    !     solve M (dxdr) = u y Pijk
    !           M (dxdz) =-u x Pijk

    !use global_parameters
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mtheta,mpsi,mgrid,pi,zetamin,deltaz,deltat,R_0
    use field_array,only: psi_grid,qtinv,Jcbn
    implicit none
#include "finclude/petscsys.h"
    type(solver_data) :: solver
    integer :: nn
    integer :: ierr
    real(wp),intent(in) :: bb_gk(nn)
    real(wp),intent(inout) :: xa_gk(nn), xt_gk(nn)
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq,mt
    real(doubleprec),allocatable,dimension(:) :: xa, xt, bb, ba, bt, fcj
    real(wp) zeta, pdum, dt, tdum, tmp
    real(wp), dimension(1) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

    real(wp),allocatable,dimension(:,:) :: xlt
    integer,dimension(:),allocatable :: ind
    integer :: fe_order, fe_ns, fe_nvar
    integer :: i,ii,j,jj,k,kk,ix,iy
    logical :: flg
    integer :: whichone = 5


    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-whichone',whichone,flg,ierr)
    if(whichone.lt.3) then
       ! in (r,z) cylindrical coordinates
       ix=1
       iy=2
    else
       ! in (a,t) flux coordinates
       ix=3
       iy=4
    endif
#ifdef CJTEST
    if(mype.eq.0) write(0, *) "	finite element: find dxda dxdt for solve",whichone,ix,iy
#endif

    fe_order=solver%fem%order
    fe_ns=3
    fe_nvar=3
    allocate( xlt(2,fe_ns) )
    allocate( ind(fe_ns) )

    neq = sum(mtheta)
    allocate(bb(neq),xa(neq),ba(neq),xt(neq),bt(neq),fcj(neq))
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi 
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             bb(ivsl) = bb_gk(ivgts)
             fcj(ivsl) = Jcbn(ivgts)/(solver%fem%points(1,ivsl)+R_0)/(solver%fem%points(1,ivsl)+R_0)
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'

    if(whichone.lt.3) fcj=1.
       !do ivsl=1,neq
       !   aaa(1)=solver%fem%points(3,ivsl)
       !   the(1)=solver%fem%points(4,ivsl)
       !   call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj, aaa(1), the(1),1)
       !   fcj(ivsl) = (zaj(1)*rqj(1)-zqj(1)*raj(1))/rj(1)
       !enddo

!   form rhs ba&bt for dxdr&dxdz respectively
    ba=0.
    bt=0.
    do ie=1, solver%fem%iemax
       ind(:)=solver%fem%iet(:,ie)
       if(ind(1) .eq. 0) stop 'ind(1) == 0'

       !if(whichone.lt.3) then
       !   xlt(1,:) = solver%fem%points(1,ind(:))
       !   xlt(2,:) = solver%fem%points(2,ind(:))
       !else
       xlt(1,:) = solver%fem%points(ix,ind(:))
       xlt(2,:) = solver%fem%points(iy,ind(:))
       if(whichone.ge.3) call fixtheta(ie,fe_order,fe_ns,xlt)
       !endif

       do i = 1,fe_nvar
          ii=solver%fem%iet(i,ie)
          do j = 1,fe_nvar
             jj=solver%fem%iet(j,ie)
          do k = 1,fe_ns
             kk=solver%fem%iet(k,ie)
             ba(ii)=ba(ii)+bb(jj)*xlt(2,k)*solver%fem%Pijk(k,j,i)*fcj(jj)
             bt(ii)=bt(ii)-bb(jj)*xlt(1,k)*solver%fem%Pijk(k,j,i)*fcj(jj)
          enddo
          enddo
       enddo
    enddo
    !call find_rhs(neq, &
    !              ba, &
    !              bt, &
    !              bb, &
    !              solver%fem%iemax, &
    !              solver%fem%iet(1,:), &
    !              solver%fem%iet(2,:), &
    !              solver%fem%iet(3,:), &
    !              solver%fem%points(1,:), &
    !              solver%fem%points(2,:))

!   solve xa&xt representing dxdr&dxdz respectively
    xa=0.
    xt=0.
    !re-enforce in case somehting went wrong
    solver%order=0
    ! 'fcj' serves as a place holder
    call petsc_solve( solver, neq, ba, xa, fcj, ierr)
    call petsc_solve( solver, neq, bt, xt, fcj, ierr)
    
    ! set solution (output)
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts+1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             xa_gk(ivgts) = xa(ivsl)
             xt_gk(ivgts) = xt(ivsl)
          else
             xa_gk(ivgts) = xa(ivsl-itmax+1)
             xt_gk(ivgts) = xt(ivsl-itmax+1)
          endif
       enddo
    enddo

    deallocate( xlt )
    deallocate( ind )
    deallocate(bb,ba,bt,xa,xt,fcj)

  end subroutine find_dxdadxdt



!========================================================================
!         find dxdr only
!========================================================================
  subroutine find_dxda( solver, nn, bb_gk, xa_gk, ierr )

    !     solve M (dxdr) = u y Pijk

    !use global_parameters
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mtheta,mpsi,mgrid,pi,zetamin,deltaz,deltat,R_0
    use field_array,only: psi_grid,qtinv,Jcbn
    implicit none
#include "finclude/petscsys.h"
    type(solver_data) :: solver
    integer :: nn
    integer :: ierr
    real(wp),intent(in) :: bb_gk(nn)
    real(wp),intent(inout) :: xa_gk(nn)
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq,mt
    real(doubleprec),allocatable,dimension(:) :: xa, xt, bb, ba, fcj
    real(wp) zeta, pdum, dt, tdum, tmp
    real(wp), dimension(1) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

    real(wp),allocatable,dimension(:,:) :: xlt
    integer,dimension(:),allocatable :: ind
    integer :: fe_order, fe_ns, fe_nvar
    integer :: i,ii,j,jj,k,kk,ix,iy
    logical :: flg
    integer :: whichone = 5


    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-whichone',whichone,flg,ierr)
    if(whichone.lt.3) then
       ! in (r,z) cylindrical coordinates
       ix=1
       iy=2
    else
       ! in (a,t) flux coordinates
       ix=3
       iy=4
    endif
#ifdef CJTEST
    if(mype.eq.0) write(0, *) "	finite element: find dxda for solve",whichone,ix,iy
#endif

    fe_order=solver%fem%order
    fe_ns=3
    fe_nvar=3
    allocate( xlt(2,fe_ns) )
    allocate( ind(fe_ns) )

    neq = sum(mtheta)
    allocate(bb(neq),xa(neq),ba(neq),fcj(neq))
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi 
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             bb(ivsl) = bb_gk(ivgts)
             fcj(ivsl) = Jcbn(ivgts)/(solver%fem%points(1,ivsl)+R_0)/(solver%fem%points(1,ivsl)+R_0)
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'

    if(whichone.lt.3) fcj=1.

!   form rhs ba for dxdr respectively
    ba=0.
    do ie=1, solver%fem%iemax
       ind(:)=solver%fem%iet(:,ie)
       if(ind(1) .eq. 0) stop 'ind(1) == 0'

       !if(whichone.lt.3) then
       !   xlt(1,:) = solver%fem%points(1,ind(:))
       !   xlt(2,:) = solver%fem%points(2,ind(:))
       !else
       xlt(1,:) = solver%fem%points(ix,ind(:))
       xlt(2,:) = solver%fem%points(iy,ind(:))
       if(whichone.ge.3) call fixtheta(ie,fe_order,fe_ns,xlt)
       !endif

       do i = 1,fe_nvar
          ii=solver%fem%iet(i,ie)
          do j = 1,fe_nvar
             jj=solver%fem%iet(j,ie)
          do k = 1,fe_ns
             kk=solver%fem%iet(k,ie)
             ba(ii)=ba(ii)+bb(jj)*xlt(2,k)*solver%fem%Pijk(k,j,i)*fcj(jj)
          enddo
          enddo
       enddo
    enddo

!   solve xa representing dxdr respectively
    xa=0.
    !re-enforce in case somehting went wrong
    solver%order=0
    ! 'fcj' serves as a place holder
    call petsc_solve( solver, neq, ba, xa, fcj, ierr)
    
    ! set solution (output)
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts+1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             xa_gk(ivgts) = xa(ivsl)
          else
             xa_gk(ivgts) = xa(ivsl-itmax+1)
          endif
       enddo
    enddo

    deallocate( xlt )
    deallocate( ind )
    deallocate(bb,ba,xa,fcj)

  end subroutine find_dxda



!========================================================================
!         find dxdz only
!========================================================================
  subroutine find_dxdt( solver, nn, bb_gk, xt_gk, ierr )

    !           M (dxdz) =-u x Pijk

    !use global_parameters
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mtheta,mpsi,mgrid,pi,zetamin,deltaz,deltat,R_0
    use field_array,only: psi_grid,qtinv,Jcbn
    implicit none
#include "finclude/petscsys.h"
    type(solver_data) :: solver
    integer :: nn
    integer :: ierr
    real(wp),intent(in) :: bb_gk(nn)
    real(wp),intent(inout) :: xt_gk(nn)
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq,mt
    real(doubleprec),allocatable,dimension(:) :: xt, bb, bt, fcj
    real(wp) zeta, pdum, dt, tdum, tmp
    real(wp), dimension(1) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

    real(wp),allocatable,dimension(:,:) :: xlt
    integer,dimension(:),allocatable :: ind
    integer :: fe_order, fe_ns, fe_nvar
    integer :: i,ii,j,jj,k,kk,ix,iy
    logical :: flg
    integer :: whichone = 5


    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-whichone',whichone,flg,ierr)
    if(whichone.lt.3) then
       ! in (r,z) cylindrical coordinates
       ix=1
       iy=2
    else
       ! in (a,t) flux coordinates
       ix=3
       iy=4
    endif
#ifdef CJTEST
    if(mype.eq.0) write(0, *) "	finite element: find dxdt for solve",whichone,ix,iy
#endif

    fe_order=solver%fem%order
    fe_ns=3
    fe_nvar=3
    allocate( xlt(2,fe_ns) )
    allocate( ind(fe_ns) )

    neq = sum(mtheta)
    allocate(bb(neq),xt(neq),bt(neq),fcj(neq))
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi 
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts + 1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             bb(ivsl) = bb_gk(ivgts)
             fcj(ivsl) = Jcbn(ivgts)/(solver%fem%points(1,ivsl)+R_0)/(solver%fem%points(1,ivsl)+R_0)
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'

    if(whichone.lt.3) fcj=1.

!   form rhs bt for dxdz respectively
    bt=0.
    do ie=1, solver%fem%iemax
       ind(:)=solver%fem%iet(:,ie)
       if(ind(1) .eq. 0) stop 'ind(1) == 0'

       !if(whichone.lt.3) then
       !   xlt(1,:) = solver%fem%points(1,ind(:))
       !   xlt(2,:) = solver%fem%points(2,ind(:))
       !else
       xlt(1,:) = solver%fem%points(ix,ind(:))
       xlt(2,:) = solver%fem%points(iy,ind(:))
       if(whichone.ge.3) call fixtheta(ie,fe_order,fe_ns,xlt)
       !endif

       do i = 1,fe_nvar
          ii=solver%fem%iet(i,ie)
          do j = 1,fe_nvar
             jj=solver%fem%iet(j,ie)
          do k = 1,fe_ns
             kk=solver%fem%iet(k,ie)
             bt(ii)=bt(ii)-bb(jj)*xlt(1,k)*solver%fem%Pijk(k,j,i)*fcj(jj)
          enddo
          enddo
       enddo
    enddo

!   solve xt representing dxdz respectively
    xt=0.
    !re-enforce in case somehting went wrong
    solver%order=0
    ! 'fcj' serves as a place holder
    call petsc_solve( solver, neq, bt, xt, fcj, ierr)
    
    ! set solution (output)
    ivsl = 0
    ivgts = 0
    do ir=0,mpsi
       itmax=mtheta(ir)
       do it=0,itmax
          ivgts = ivgts+1
          if(it.lt.itmax) then ! skip double stored 0,2pi value
             ivsl = ivsl + 1
             xt_gk(ivgts) = xt(ivsl)
          else
             xt_gk(ivgts) = xt(ivsl-itmax+1)
          endif
       enddo
    enddo

    deallocate(bb,bt,xt,fcj)
    deallocate( xlt )
    deallocate( ind )

  end subroutine find_dxdt



!========================================================================
!         fix theta periodic behavior
!========================================================================
  subroutine fixtheta(ie,fe_order,fe_ns,xlt)
    use global_parameters, only:mype,mgrid,pi
    implicit none
    integer :: ie,fe_order,fe_ns
    real(wp),intent(inout) :: xlt(2,fe_ns)
    real(wp) :: eps, big

    eps=pi/6.
    big=pi*1.5

    if(fe_order.eq.2) then
        if( xlt(2,1).lt.eps .and. xlt(2,3).gt.big .and. xlt(2,5).gt.big ) then
            xlt(2,1) = 2.*pi + xlt(2,1)
            !if(mype .eq. 0) write(0,*) '	t1 changed:', ie, xlt(2,1)
        endif
        if( xlt(2,3).lt.eps .and. xlt(2,5).gt.big .and. xlt(2,1).gt.big ) then
            xlt(2,3) = 2.*pi + xlt(2,3)
            !if(mype .eq. 0) write(0,*) '	t2 changed:', ie, xlt(2,3)
        endif
        if( xlt(2,5).lt.eps .and. xlt(2,1).gt.big .and. xlt(2,3).gt.big ) then
            xlt(2,5) = 2.*pi + xlt(2,5)
            !if(mype .eq. 0) write(0,*) '	t3 changed:', ie, xlt(2,5)
        endif
        if( xlt(2,1).lt.eps .and. xlt(2,3).lt.eps .and. xlt(2,5).gt.big ) then
            xlt(2,1) = 2.*pi + xlt(2,1)
            xlt(2,3) = 2.*pi + xlt(2,3)
            !if(mype .eq. 0) write(0,*) '	t4 changed:', ie, xlt(2,1), xlt(2,3)
        endif
        if( xlt(2,3).lt.eps .and. xlt(2,5).lt.eps .and. xlt(2,1).gt.big ) then
            xlt(2,3) = 2.*pi + xlt(2,3)
            xlt(2,5) = 2.*pi + xlt(2,5)
            !if(mype .eq. 0) write(0,*) '	t5 changed:', ie, xlt(2,3), xlt(2,5)
        endif
        if( xlt(2,5).lt.eps .and. xlt(2,1).lt.eps .and. xlt(2,3).gt.big ) then
            xlt(2,5) = 2.*pi + xlt(2,5)
            xlt(2,1) = 2.*pi + xlt(2,1)
            !if(mype .eq. 0) write(0,*) '	t6 changed:', ie, xlt(2,5), xlt(2,1)
        endif
        xlt(:,2) = .5*(xlt(:,1)+xlt(:,3))
        xlt(:,4) = .5*(xlt(:,3)+xlt(:,5))
        xlt(:,6) = .5*(xlt(:,5)+xlt(:,1))
    else
        if( xlt(2,1).lt.eps .and. xlt(2,2).gt.big .and. xlt(2,3).gt.big ) then
            xlt(2,1) = 2.*pi + xlt(2,1)
            !if(mype .eq. 0) write(0,*) '	t1 changed:', ie, xlt(2,1)
        endif
        if( xlt(2,2).lt.eps .and. xlt(2,3).gt.big .and. xlt(2,1).gt.big ) then
            xlt(2,2) = 2.*pi + xlt(2,2)
            !if(mype .eq. 0) write(0,*) '	t2 changed:', ie, xlt(2,2)
        endif
        if( xlt(2,3).lt.eps .and. xlt(2,1).gt.big .and. xlt(2,2).gt.big ) then
            xlt(2,3) = 2.*pi + xlt(2,3)
            !if(mype .eq. 0) write(0,*) '	t3 changed:', ie, xlt(2,3)
        endif
        if( xlt(2,1).lt.eps .and. xlt(2,2).lt.eps .and. xlt(2,3).gt.big ) then
            xlt(2,1) = 2.*pi + xlt(2,1)
            xlt(2,2) = 2.*pi + xlt(2,2)
            !if(mype .eq. 0) write(0,*) '	t4 changed:', ie, xlt(2,1), xlt(2,2)
        endif
        if( xlt(2,2).lt.eps .and. xlt(2,3).lt.eps .and. xlt(2,1).gt.big ) then
            xlt(2,2) = 2.*pi + xlt(2,2)
            xlt(2,3) = 2.*pi + xlt(2,3)
            !if(mype .eq. 0) write(0,*) '	t5 changed:', ie, xlt(2,2), xlt(2,3)
        endif
        if( xlt(2,3).lt.eps .and. xlt(2,1).lt.eps .and. xlt(2,2).gt.big ) then
            xlt(2,3) = 2.*pi + xlt(2,3)
            xlt(2,1) = 2.*pi + xlt(2,1)
            !if(mype .eq. 0) write(0,*) '	t6 changed:', ie, xlt(2,3), xlt(2,1)
        endif
    endif
  end subroutine fixtheta


#ifdef OLDWAY
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
    !   call VecCreateSeqWithArray(MPI_COMM_SELF,nn,yy,vec_yy,ierr)
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

#else
!========================================================================
!         petsc_solve
!========================================================================
  subroutine petsc_solve( solver, nn, bb, xx, bh, ierr )
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
    real(doubleprec),intent(in) :: bb(nn), bh(nn)
    real(doubleprec),intent(inout) :: xx(nn)
    type(solver_data) :: solver
    ! locals
    Vec       vec_yy,rhsVec,tmpVec01,tmpVec10,vec_xx,vec_bh,tmpVecbh
    integer          iplot,i,j
    PetscScalar t1
    PetscScalar,allocatable,dimension(:) :: yy, zz
    PetscViewer viewer
    integer :: m,n,p,vsize,msize,nsize
    PetscInt :: its
    PetscReal :: rnorm


    if(solver%mat%data .eq. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    !     copy real vector into double precision
    allocate(yy(nn),zz(nn))

    !     scatter rhs into PETSc vector
    yy = bb
    call VecCreateSeqWithArray(MPI_COMM_SELF,1,nn,yy,vec_yy,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    !     call VecSet( solver%bVec, 1.e30, ierr ) ! debug
    call VecScatterBegin(solver%scat,vec_yy, solver%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    call VecScatterEnd(solver%scat,vec_yy, solver%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
    if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

    if(solver%order .ge. 1) then 
       !cj M * b =>> rhs for solve Stiff xx = rhs
       call VecDuplicate(solver%bVec, rhsVec, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
       call VecCopy(solver%bVec, rhsVec, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
       call MatMult(solver%mat%mdata, rhsVec, solver%bVec, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
       !  call VecView( solver%bVec, PETSC_VIEWER_STDOUT_WORLD, ierr )

       !! add bh and neumann bc fix to M*bb to M*bb
       !  if(mype==0) write(*,*)"	petsc_solve: add hh to bb"
       call VecDuplicate(solver%bVec, tmpVecbh, ierr)
       call VecCreateSeqWithArray(MPI_COMM_SELF,1,nn,bh,vec_bh,ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
       call VecScatterBegin(solver%scat,vec_bh, tmpVecbh, INSERT_VALUES,SCATTER_FORWARD,ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
       call VecScatterEnd(solver%scat,vec_bh, tmpVecbh, INSERT_VALUES,SCATTER_FORWARD,ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
       call VecAXPBY(solver%bVec, 1., 1., tmpVecbh, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
       call VecDestroy( vec_bh, ierr )
       call VecDestroy( tmpVecbh, ierr )

       ! apply dirichlet bc
       ! bVec = bVec * zero10 + xVec * zero01
       call VecCreateSeqWithArray(MPI_COMM_SELF,1,nn,xx,vec_xx,ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
       call VecScatterBegin(solver%scat,vec_xx, solver%xVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
       call VecScatterEnd(solver%scat,vec_xx, solver%xVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

       ! take the part of bVec corresponding to internal vertex
       call VecDuplicate(solver%bVec, tmpVec10, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
       call VecPointwiseMult(tmpVec10, solver%bVec, solver%zeroVec10, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

       ! take the part of xVec corresponding to boundary vertex
       call VecDuplicate(solver%bVec, tmpVec01, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
       call VecPointwiseMult(tmpVec01, solver%xVec, solver%zeroVec01, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )

       call VecAXPBYPCZ(solver%bVec, 1., 1., 0., tmpVec10, tmpVec01, ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr )
    endif
    !   call VecGetSize( solver%bVec, vsize, ierr )
    !   if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    !       if(mype==0) write(*,*) 'dimension=',nn, 'vsize=',vsize
    !   call MatGetSize(solver%mat%mdata,msize,nsize,ierr)
    !   if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    !       if(mype==0) write(*,*) 'mdata: msize=',msize, 'nsize=',nsize
    !   call MatGetSize(solver%mat%data,msize,nsize,ierr)
    !   if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    !       if(mype==0) write(*,*) 'data: msize=',msize, 'nsize=',nsize

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
    call VecSet(solver%xVec, 0., ierr)
    call VecSet(solver%bVec, 0., ierr)
    if(solver%order .ge. 1) then
       call VecDestroy( vec_xx, ierr )
       call VecDestroy( rhsVec, ierr )
       call VecDestroy( tmpVec01, ierr )
       call VecDestroy( tmpVec10, ierr )
    endif
  end subroutine petsc_solve
#endif


!========================================================================
!         set up gk_matid for Dirichlet boundary condition
!========================================================================
  subroutine new_gk_matid( this, bc_type )
    use precision,only:wp,doubleprec
    use global_parameters,only: mype,mtheta,mpsi,ZF_SEPARATION,ifelectron
    use particle_decomp,only:myrank_partd,nproc_partd
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
    type(phit_matrix_data) :: this
    integer bc_type
    ! locals
    integer i, j, ierr, it, id

    ! form map from gts node numbering (includes 0 & 2pi) to matrix numbering
    it = sum(mtheta) !... + mpsi + 2
    allocate(this%gk_matid(it))
    this%gk_matid = 0
    !... this%gk_matid(it) = -2
    
    id = 0
    it = 0
    do i=0,mpsi
       do j=0,mtheta(i)
          it = it + 1
          if( j == mtheta(i) ) then
!            this%gk_matid(it) = id - mtheta(i) + 1 - 1
          else
             id = id + 1
             this%gk_matid(id) = id

             ! BCs
             if(bc_type==1) then      ! for Poisson solver
                ! Dirichlet boundary condiction at inner boundary: phi(i=0)=0
                if( i==0 ) this%gk_matid(id) = -1
             else if(bc_type==2) then      ! for dphi/dt solver
                if( i==0 ) this%gk_matid(id) = -1
                if( i==mpsi ) this%gk_matid(id) = -1
             else if(bc_type==3) then
                ! neumann bc at a0, dirichlet bc at a1
                if( i==mpsi ) this%gk_matid(id) = -1
             else                      ! for 
                ! bc_type = 0
                ! neumann bc at both a0 and a1
             endif
             !...if(mype .eq. 0 ) write(112,112) id, i, j, this%gk_matid(id)
             !...112 format("new_gk_matid ", 4(2x,i7))
          end if
       end do
    end do
    if(id .ne. sum(mtheta)) call MPI_ABORT(mpi_comm_world,1,ierr)
    !... if(this%gk_matid(sum(mtheta) + mpsi + 2) .ne. -2) call MPI_ABORT(mpi_comm_world,1,ierr)
  end subroutine new_gk_matid


!========================================================================
!         rhs fix for Neumann boundary condition
!========================================================================
subroutine neumbcfix(i,mtheta_i,offset,neq,fem,dphida,dphidt,bn)
    use precision,only:wp,doubleprec
    use fem_data_class
    use esi_interface, only : metric
    implicit none
#include "finclude/petscsys.h"
#include "finclude/petscmat.h" 
#include "finclude/petscvec.h" 
#include "finclude/petscviewer.h"
    integer i,mtheta_i,offset,neq
    type(fem_data) :: fem
    real(wp) :: dphida,dphidt
    real(wp) :: bn(neq)
    ! locals
    real(wp) :: gs,dphida1,dphidt1,dphida2,dphidt2,x1,x2,y1,y2,l12
    real(wp) :: gaa1, gtt1, gat1, gaa2, gtt2, gat2
    real(wp), dimension(2) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj, fcj
    integer :: j,i1,i2,mt
    integer, external :: metric_elements

       !set gs, dphida, dphidt as a place holder for input into "solve" subroutine
       gs=1.
       do j=1,mtheta_i
          dphida1=1.  !to be replaced by dphida at j
          dphidt1=1.  !to be replaced by dphidt at j
          dphida2=1.  !to be replaced by dphida at j
          dphidt2=1.  !to be replaced by dphidt at j

          ! the boundary edge is make of the following 2 vertices
          i1 = offset + j
          i2 = offset + j+1
          ! if j is the last index, then i1 connects to the first index
          if(j.eq.mtheta_i) i2=offset + 1

          aaa(1)=fem%points(3,i1)
          the(1)=fem%points(4,i1)
          aaa(2)=fem%points(3,i2)
          the(2)=fem%points(4,i2)
          call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,2)
          fcj = zaj*rqj-zqj*raj/rj

          ! if (it, ir) makes one edge
          x1=fem%points(3,i1)
          y1=fem%points(4,i1)
          x2=fem%points(3,i2)
          y2=fem%points(4,i2)
          l12=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)

          mt=metric_elements(aaa(1),the(1),0.)
          gaa1=metric(1)
          gtt1=metric(2)
          gat1=metric(4)
          bn(i1) = bn(i1) - .5*fcj(1)*gs*(gaa1*dphida1)*sqrt(l12)
          bn(i1) = bn(i1) - .5*fcj(1)*gs*(gat1*dphidt1)*sqrt(l12)

          mt=metric_elements(aaa(2),the(2),0.)
          gaa2=metric(1)
          gtt2=metric(2)
          gat2=metric(4)
          bn(i2) = bn(i2) - .5*fcj(2)*gs*(gaa2*dphida2)*sqrt(l12)
          bn(i2) = bn(i2) - .5*fcj(2)*gs*(gat2*dphidt2)*sqrt(l12)
       end do
end subroutine neumbcfix



!========================================================================
!         add_hh2bb
!========================================================================
  subroutine add_hh2bb( fem, nn, hh, bh, ierr )
    !     add hh(dphi/dt) term to bb(rhs) before KSPSolve
    !     
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mtheta,mpsi,mgrid,zetamin,deltaz,deltat,R_0,a0,a1,pi,csc0,m0ms,c0cs_i,msm0_i
    use field_array,only: psi_grid,qtinv,Jcbn,gij,gsij
    use esi_interface, only : metric
    use fem_data_class
    implicit none
#include "finclude/petscsys.h"
    type(fem_data) :: fem
    integer,intent(in) :: nn
    real(wp),intent(in) :: hh(nn)
    real(wp),intent(inout) :: bh(nn) 
    integer,intent(out) :: ierr
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq
    real(wp), dimension(1) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
    real(wp) :: fcj,beta,gs
    real(wp) zeta, pdum, dt, tdum
    integer :: mt,i, ii, j, jj, k, kk
    real(doubleprec),allocatable,dimension(:) :: bgsjgaa,bgsjgat,bgsjgaa_ho,bgsjgat_ho
    real(wp) alpha, beta, ns0, ne0, gs, Zs, gaa, gat
    real(wp), external ::temperature_e, density
    integer, external :: metric_elements
    integer :: i0, i1, i2, i3, i4, i5
    integer :: m,mv,n,nv
    real(wp),allocatable,dimension(:,:) :: xlt
    integer,dimension(:),allocatable :: ind
    integer,allocatable,dimension(:) :: indlocal, indglobal
    integer :: fe_order, fe_ns, fe_nvar, fe_ndof
    logical :: flg
    integer :: whichone = 5


   !if(mype==0) write(*,*)"	finite element: add hh to bb fem_order=", fem%order
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-whichone',whichone,flg,ierr)

    fe_order=fem%order
    fe_ns=3
    fe_ndof=1
    if(fe_order .gt. 2) fe_ndof=6
    if(fe_order .eq. 2) fe_ns=6
    fe_nvar=fe_ns*fe_ndof
    allocate( xlt(2,fe_ns) )
    allocate( ind(fe_ns) )
    allocate(indlocal(fe_nvar), indglobal(fe_nvar))

    !cal beta * gs * fluxcoordinatejacobian * gradadotgrada or gradadotgradt
    neq = sum(mtheta)
    allocate(bgsjgaa(neq*fe_ndof), bgsjgat(neq*fe_ndof))
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
             tdum=tdum/(2.0_wp*pi)
             tdum=2.0_wp*pi*(tdum-aint(tdum))

             aaa(1)=fem%points(3,ivsl)
             the(1)=fem%points(4,ivsl)
             ns0=density(aaa(1),0.0_wp)
             ne0=density(aaa(1),0.0_wp)
             beta = - 1./ ne0

             gs = gsij(ivgts)
             gaa = gij(1,ivgts)
             gat = gij(4,ivgts)
             fcj = Jcbn(ivgts)/(fem%points(1,ivsl)+R_0)/(fem%points(1,ivsl)+R_0)

             i0=(ivsl-1)*fe_ndof+1
             bgsjgaa(i0) = beta*gs*fcj*gaa
             bgsjgat(i0) = beta*gs*fcj*gat
             !if(mype.eq.0)write(117,117) ivsl, gs, gaa, gat, bgsjgaa(i0), bgsjgat(i0)
             !117 format("addhh2bb=", 1(2x,i5), 5(2x,d12.6))
             do k=2,fe_ndof
                i0=(ivsl-1)*fe_ndof+k
                bgsjgaa(i0) = beta*gs*fcj*gaa
                bgsjgat(i0) = beta*gs*fcj*gat
             enddo
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    
    if( fem%order .eq. 2) then
        ! interpolate linear to ho array
        allocate(bgsjgaa_ho(fem%nmesh_pts), bgsjgat_ho(fem%nmesh_pts))
        call hogrid(fem%iemax, fem%iet, neq, bgsjgaa, fem%nmesh_pts, bgsjgaa_ho)
        call hogrid(fem%iemax, fem%iet, neq, bgsjgat, fem%nmesh_pts, bgsjgat_ho)

        !cal extra term to be added to M*bb(rhs)
        do ie=1, fem%iemax
           do i = 1,fe_ns
              ii=fem%iet(i,ie)
              do j = 1,fe_ns
                 jj=fem%iet(j,ie)
              do k = 1,fe_ns
                 kk=fem%iet(k,ie)
                 bh(ii)=bh(ii)+hh(jj)*fem%points(4,kk)*fem%Pijk(k,i,j)*bgsjgaa_ho(jj)
                 bh(ii)=bh(ii)-hh(jj)*fem%points(3,kk)*fem%Pijk(k,i,j)*bgsjgat_ho(jj)
              enddo
              enddo
           enddo
        enddo

        deallocate(bgsjgaa_ho, bgsjgat_ho)

    else
        !cal extra term to be added to M*bb(rhs)
        do ie=1, fem%iemax
           ind(:)=fem%iet(:,ie)
           if(ind(1) .eq. 0) stop 'ind(1) == 0'

           if(whichone.lt.3) then
              ! in (r,z) cylindrical coordinates
              xlt(1,:) = fem%points(1,ind(:))
              xlt(2,:) = fem%points(2,ind(:))
           else
              ! in (a,t) flux coordinates
              xlt(1,:) = fem%points(3,ind(:))
              xlt(2,:) = fem%points(4,ind(:))
              call fixtheta(ie,fe_order,fe_ns,xlt)
           endif

           !preparing for index
           do m = 1,fe_ns
              mv=fem%iet(m,ie)
              do n = 1,fe_ndof
                 nv=(m-1)*fe_ndof+n
                 indlocal(nv)=nv
                 indglobal(nv)=(mv-1)*fe_ndof+n
              enddo
           enddo

           do i = 1,fe_nvar
              ii=indglobal(i)
              do j = 1,fe_nvar
                 jj=indglobal(j)
              do k = 1,fe_ns
                 kk=fem%iet(k,ie)
                 bh(ii)=bh(ii)-hh(jj)*xlt(2,k)*fem%Pijk(k,j,i)*bgsjgaa(jj)
                 bh(ii)=bh(ii)+hh(jj)*xlt(1,k)*fem%Pijk(k,j,i)*bgsjgat(jj)
              enddo
              enddo
           enddo
        enddo
    endif

    deallocate( xlt )
    deallocate( ind )
    deallocate(indlocal, indglobal)
    deallocate(bgsjgaa, bgsjgat)
  end subroutine add_hh2bb



!========================================================================
!         add_hh2bb
!========================================================================
  subroutine add_hh2bb_0( fem, nn, hh, bh, ierr )
    !     add hh(dphi/dt) term to bb(rhs) before KSPSolve
    !     
    use precision,only:wp,doubleprec
    use global_parameters, only:mype,mtheta,mpsi,mgrid,zetamin,deltaz,deltat,R_0,a0,a1,pi,csc0,m0ms,c0cs_i,msm0_i
    use field_array,only: psi_grid,qtinv,Jcbn,gij,gsij
    use esi_interface, only : metric
    use fem_data_class
    implicit none
#include "finclude/petscsys.h"
    type(fem_data) :: fem
    integer,intent(in) :: nn
    real(wp),intent(in) :: hh(nn)
    real(wp),intent(inout) :: bh(nn) 
    integer,intent(out) :: ierr
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq
    real(wp), dimension(1) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
    real(wp) :: fcj,beta,gs
    real(wp) zeta, pdum, dt, tdum
    integer :: mt,i, ii, j, jj, k, kk
    real(doubleprec),allocatable,dimension(:) :: bgsjgaa,bgsjgat,bgsjgaa_ho,bgsjgat_ho
    real(wp) alpha, beta, ns0, ne0, gs, Zs, gaa, gat
    real(wp), external ::temperature_e, density
    integer, external :: metric_elements
    integer :: i0, i1, i2, i3, i4, i5
    integer :: m,mv,n,nv
    real(wp),allocatable,dimension(:,:) :: xlt
    integer,dimension(:),allocatable :: ind
    integer,allocatable,dimension(:) :: indlocal, indglobal
    integer :: fe_order, fe_ns, fe_nvar, fe_ndof
    logical :: flg
    integer :: whichone = 5


   !if(mype==0) write(*,*)"	finite element: add hh to bb fem_order=", fem%order
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-whichone',whichone,flg,ierr)

    fe_order=fem%order
    fe_ns=3
    fe_ndof=1
    if(fe_order .gt. 2) fe_ndof=6
    if(fe_order .eq. 2) fe_ns=6
    fe_nvar=fe_ns*fe_ndof
    allocate( xlt(2,fe_ns) )
    allocate( ind(fe_ns) )
    allocate(indlocal(fe_nvar), indglobal(fe_nvar))

    !cal beta * gs * fluxcoordinatejacobian * gradadotgrada or gradadotgradt
    neq = sum(mtheta)
    allocate(bgsjgaa(neq*fe_ndof), bgsjgat(neq*fe_ndof))
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
             tdum=tdum/(2.0_wp*pi)
             tdum=2.0_wp*pi*(tdum-aint(tdum))

             aaa(1)=fem%points(3,ivsl)
             the(1)=fem%points(4,ivsl)
             ns0=density(aaa(1),0.0_wp)
             ne0=density(aaa(1),0.0_wp)
             beta = - 1./ ne0

             gs = gsij(ivgts)
             gaa = gij(1,ivgts)
             gat = gij(4,ivgts)
             fcj = Jcbn(ivgts)/(fem%points(1,ivsl)+R_0)/(fem%points(1,ivsl)+R_0)

             i0=(ivsl-1)*fe_ndof+1
             bgsjgaa(i0) = beta*gs*fcj*gaa
             bgsjgat(i0) = beta*gs*fcj*gat
             !if(mype.eq.0)write(117,117) ivsl, gs, gaa, gat, bgsjgaa(i0), bgsjgat(i0)
             !117 format("addhh2bb=", 1(2x,i5), 5(2x,d12.6))
             do k=2,fe_ndof
                i0=(ivsl-1)*fe_ndof+k
                bgsjgaa(i0) = beta*gs*fcj*gaa
                bgsjgat(i0) = beta*gs*fcj*gat
             enddo
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    
    if( fem%order .eq. 2) then
        ! interpolate linear to ho array
        allocate(bgsjgaa_ho(fem%nmesh_pts), bgsjgat_ho(fem%nmesh_pts))
        call hogrid(fem%iemax, fem%iet, neq, bgsjgaa, fem%nmesh_pts, bgsjgaa_ho)
        call hogrid(fem%iemax, fem%iet, neq, bgsjgat, fem%nmesh_pts, bgsjgat_ho)

        !cal extra term to be added to M*bb(rhs)
        do ie=1, fem%iemax
           do i = 1,fe_ns
              ii=fem%iet(i,ie)
              do j = 1,fe_ns
                 jj=fem%iet(j,ie)
              do k = 1,fe_ns
                 kk=fem%iet(k,ie)
                 bh(ii)=bh(ii)+hh(jj)*fem%points(4,kk)*fem%Pijk(k,i,j)*bgsjgaa_ho(jj)
                 bh(ii)=bh(ii)-hh(jj)*fem%points(3,kk)*fem%Pijk(k,i,j)*bgsjgat_ho(jj)
              enddo
              enddo
           enddo
        enddo

        deallocate(bgsjgaa_ho, bgsjgat_ho)

    else
        !cal extra term to be added to M*bb(rhs)
        do ie=1, fem%iemax
           ind(:)=fem%iet(:,ie)
           if(ind(1) .eq. 0) stop 'ind(1) == 0'

           if(whichone.lt.3) then
              ! in (r,z) cylindrical coordinates
              xlt(1,:) = fem%points(1,ind(:))
              xlt(2,:) = fem%points(2,ind(:))
           else
              ! in (a,t) flux coordinates
              xlt(1,:) = fem%points(3,ind(:))
              xlt(2,:) = fem%points(4,ind(:))
              call fixtheta(ie,fe_order,fe_ns,xlt)
           endif

           !preparing for index
           do m = 1,fe_ns
              mv=fem%iet(m,ie)
              do n = 1,fe_ndof
                 nv=(m-1)*fe_ndof+n
                 indlocal(nv)=nv
                 indglobal(nv)=(mv-1)*fe_ndof+n
              enddo
           enddo

           do i = 1,fe_nvar
              ii=indglobal(i)
              do j = 1,fe_nvar
                 jj=indglobal(j)
              do k = 1,fe_ns
                 kk=fem%iet(k,ie)
                 bh(ii)=bh(ii)-hh(jj)*xlt(2,k)*fem%Pijk(k,j,i)*bgsjgaa(jj)
                 bh(ii)=bh(ii)+hh(jj)*xlt(1,k)*fem%Pijk(k,j,i)*bgsjgat(jj)
              enddo
              enddo
           enddo
        enddo
    endif

    deallocate( xlt )
    deallocate( ind )
    deallocate(indlocal, indglobal)
    deallocate(bgsjgaa, bgsjgat)
  end subroutine add_hh2bb_0


!========================================================================
!                         subroutine form_tangent
!========================================================================
!subroutine form_tangent( Amat, Mmat, iemax, iet1, iet2, iet3, npts, points, &
!subroutine form_tangent( order, Amat, Mmat, iemax, iet, npts, points, &
subroutine form_tangent( fem, Amat, Mmat, npts, &
     alpha, beta, ierr )
  !     
  !     input:
  !     iemax: number of elements
  !     iet1,2,3[iemax]: nodes of triangles of mesh
  !     npts: number of mesh points
  !     points[6,npts]: coordinates of points in cartesian and flux coordinate
  !     
  !     output:
  !     Amat: stiffness matrix
  !     Mmat: mass matrix
  !     ierr: error code
  !
  use precision,only:wp,doubleprec
  use particle_decomp,only:myrank_partd,nproc_partd
  use global_parameters, only:mype, csc0, m0ms, c0cs_i, msm0_i, ifelectron, pi,ifdphi
  use esi_interface, only : metric
  use fem_data_class
  implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
  type(fem_data) :: fem
  integer order,iemax,npts,ierr
! integer,dimension(:),allocatable :: iet1,iet2,iet3
! integer,dimension(:,:) :: iet
! real(wp) :: points(6,npts)
  real(wp) alpha, beta
  Mat  Amat, Mmat
!  real(doubleprec)  xlt(2,3),sarr(3,3),sarr2(3*3),ul(3,4),eps ! therm2d uses *8
!  real(doubleprec)           marr(3,3),marr2(3*3)             ! mass matrix
!  real(doubleprec)  rrarr(3,3),rzarr(3,3),zrarr(3,3),zzarr(3,3),pparr(3,3),rarr(3,3),zarr(3,3)
  real(doubleprec) eps,big ! therm2d uses *8
  real(doubleprec),dimension(:,:,:),allocatable ::  fe_Pijk
  real(doubleprec),dimension(:,:),allocatable ::  xlt,sarr,marr,ul ! therm2d uses *8
  real(doubleprec),dimension(:),allocatable ::  sarr2,marr2            ! mass matrix
  real(doubleprec),dimension(:,:),allocatable ::  rrarr,rzarr,zrarr,zzarr,pparr,rarr,zarr
  real(doubleprec),dimension(:),allocatable ::  rarr2,zarr2            ! r&z matrix
  real(doubleprec) :: c1=1., c2=1.
  real(doubleprec) :: fcj=1.  ! flux coordinate jacobian
  real(doubleprec) :: cf1=1., cf2=0., cf3=1., cf4=0., cf5=0., cf6=0.
! integer ie,i,j,k,kk,nz,ind2(3),ind3(3),ind(3),nel,nst
  integer ie,i,j,k,n,kk,ki,nz,nel,nst
  integer,dimension(:),allocatable :: ind2,ind3,ind
! real(wp), dimension(3) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
! real(wp), dimension(3) :: gaa,gtt,gpp,gat,vthe3,q,ns0,ne0
  real(wp), dimension(:),allocatable :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
  real(wp), dimension(:),allocatable :: gaa,gtt,gpp,gat,vthe3,q,ns0,ne0
  real(wp) :: Zs, gs,qq

  real(wp), external ::temperature_e, density
  integer, external :: metric_elements

  ! Lagrange finite elements
  integer mt
  integer fe_order, fe_node, fe_ns, fe_nw, fe_ndof, fe_nvar
  real(doubleprec),dimension(:),allocatable :: fe_weight
  real(doubleprec),dimension(:,:),allocatable :: fe_NN
  real(doubleprec),dimension(:,:,:),allocatable :: fe_dN1, fe_Eij, fe_bEij
  logical :: flg
  integer :: whichone = 5


  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-whichone',whichone,flg,ierr)

  !cj2014.05.02 c1 = alpha
  !cj2014.05.02 c2 = beta
  eps = pi/6.
  big = pi*1.5

  call MatZeroEntries(Amat,ierr)
  call MatZeroEntries(Mmat,ierr)

  fe_order=fem%order
  call lagrange_element_define(fe_order, &
                               fe_node, &
                               fe_ns, &
                               fe_nw, &
                               fe_ndof, &
                               fe_nvar)

  nel = fe_ns                   ! points per element
  nst = fe_ns                   ! size of elem matrix
  allocate(fe_weight(fe_nw))
  allocate(fe_NN(fe_nvar,fe_nw))
  allocate(fe_dN1(fe_node,fe_nvar,fe_nw))
  allocate(fe_Eij(fe_nvar,fe_ns,fe_nw))
  allocate(fem%Pijk(fe_ns,fe_nvar,fe_nvar))
  allocate(fe_Pijk(fe_ns,fe_nvar,fe_nvar))
  allocate(fe_bEij(fe_ns,fe_ns,fe_nw))
  call lagrange_element_construct(fe_order, &
                                   fe_ns, &
                                   fe_nw, &
                                   fe_ndof, &
                                   fe_nvar, &
                                   fe_weight, &
                                   fe_NN, &
                                   fe_dN1, &
                                   fe_Eij, &
                                   fem%Pijk, &
                                   fe_bEij)

  allocate( ind2(fe_ns),ind3(fe_ns),ind(fe_ns) )
  allocate( xlt(2,fe_ns) )
  allocate( sarr(fe_nvar,fe_nvar),sarr2(fe_nvar*fe_nvar),ul(fe_nvar,4) )
  allocate( marr(fe_nvar,fe_nvar),marr2(fe_nvar*fe_nvar) )
  allocate( rrarr(fe_nvar,fe_nvar),rzarr(fe_nvar,fe_nvar),zrarr(fe_nvar,fe_nvar),zzarr(fe_nvar,fe_nvar),pparr(fe_nvar,fe_nvar) )
  allocate( rarr(fe_nvar,fe_nvar),rarr2(fe_nvar*fe_nvar) )
  allocate( zarr(fe_nvar,fe_nvar),zarr2(fe_nvar*fe_nvar) )
  ul = 0.                   ! displacement, velocit and acceloration (not used)

  allocate( aaa(fe_ns),the(fe_ns),Fj(fe_ns),Faj(fe_ns),gFaj(fe_ns),gYaj(fe_ns),rj(fe_ns),raj(fe_ns),rqj(fe_ns),zj(fe_ns),zaj(fe_ns),zqj(fe_ns),Bj(fe_ns) )
  allocate( gaa(fe_ns),gtt(fe_ns),gpp(fe_ns),gat(fe_ns),vthe3(fe_ns),q(fe_ns),ns0(fe_ns),ne0(fe_ns) )

  !july 31, 2013 to calculate dxdr, dxdz
  !shouldn't do MPI_ALLREDUCE, but just for now to fix bug
  fe_Pijk(:,:,:)= fem%Pijk(:,:,:)
    call MPI_ALLREDUCE(fe_Pijk,fem%Pijk,fe_ns*fe_nvar*fe_nvar,MPI_DOUBLE_PRECISION,MPI_SUM,partd_comm,ierr)
#ifdef CJTEST
  do i=1,fe_nvar
     do j=1,fe_nvar
        if(mype==0) write(*,115) i,j,fem%Pijk(1,j,i),fem%Pijk(2,j,i),fem%Pijk(3,j,i)
        115 format("	Pijk:", 2(2x,i4), 3(1x,d12.6))
     enddo
  enddo
#endif

  if(mype .eq. 0) write(0,*) '	form_tangent: iemax=', fem%iemax
  do ie=1,fem%iemax
     ind(:)=fem%iet(:,ie)
     if(ind(1) .eq. 0) stop 'ind(1) == 0'
    
     if(whichone.lt.3) then
        ! in (r,z) cylindrical coordinates
        xlt(1,:) = fem%points(1,ind(:))
        xlt(2,:) = fem%points(2,ind(:))
        fcj=1.
     else
        ! in (a,t) flux coordinates
        xlt(1,:) = fem%points(3,ind(:))
        xlt(2,:) = fem%points(4,ind(:))
        call fixtheta(ie,fe_order,fe_ns,xlt)
        call cmp_phy_qty
     endif

        sarr = 0.; marr = 0.
        rrarr = 0.; rzarr = 0.; zrarr = 0.; zzarr = 0.; pparr = 0.
        rarr = 0.; zarr = 0.
#ifdef USEGTSFE
        if(mype==0 .and. ie==1 ) &
        write(0,*) "	form_tangent: use historical gts elements whichone,fcj=", whichone, fcj
        call therm2d(c1,c2,ul,xlt,sarr,marr,rrarr,rzarr,zrarr,zzarr,rarr,zarr,nel,nst,2,1)
#else
        if(mype==0 .and. ie==1 ) &
        write(0,*) "	form_tangent: use lagrange elements whichone,fcj=", whichone, fcj
        call lagrange_element_matrix(fe_ns, fe_nw, fe_ndof, fe_nvar, &
                                     xlt(1,:), xlt(2,:),  &
                                     fe_weight, fe_NN, fe_dN1, fe_Eij, fe_bEij, &
                                     marr,sarr,rrarr,rzarr,zrarr,zzarr,pparr,rarr,zarr)
#endif

#ifdef NOC1ELE
        nz = 0
        do i=1,nel
           k = ind(i)
           nz = nz + 1
!          ind3(nz) = i
           ind2(nz) = k - 1 ! zero based
        enddo
#else
        nz = 0
        do i=1,nel
        do ki=1,fe_ndof
           k = ind(i)*fe_ndof + ki
           nz = nz + 1
           ind2(nz) = k - 1*fe_ndof - 1 ! zero based
        enddo
        enddo
#endif

        kk = 0
        do j=1,nz
           k = ind(i)
           do i=1,nz
             kk = kk + 1
                if(whichone==1) then
             sarr2(kk) = cf1*rrarr(i,j) &
                        +cf2*2.*rzarr(i,j) &
                        +cf3*zzarr(i,j) &
                        +cf4*rarr(i,j) &
                        +cf5*zarr(i,j)
                else if(whichone==2) then
             sarr2(kk) = sarr(i,j)
                else if(whichone==3) then
             sarr2(kk) = fcj*( rrarr(i,j) + 2.*rzarr(i,j) + zzarr(i,j) )
                else if(whichone==4) then
             sarr2(kk) = fcj*( cf1*rrarr(i,j) + 2.*cf2*rzarr(i,j) + cf3*zzarr(i,j) )
                else
             sarr2(kk) = alpha * fcj * marr(i,j) &
                       + beta * fcj * gs * ( cf1*rrarr(i,j) &
                       +                     cf2*rzarr(i,j) &
                       +                     cf2*zrarr(i,j) &
                       +                     cf3*zzarr(i,j) &
                       +                     cf4*pparr(i,j)/qq &
                                           )
                endif
             marr2(kk) = fcj*marr(i,j)
            !rarr2(kk) = fcj*rarr(i,j)
            !zarr2(kk) = fcj*zarr(i,j)
           enddo
        enddo

        call MatSetValues(Amat,nz,ind2,nz,ind2,sarr2,ADD_VALUES,ierr)
        call MatSetValues(Mmat,nz,ind2,nz,ind2,marr2,ADD_VALUES,ierr)
  enddo

  deallocate( ind2,ind3,ind )
  deallocate( xlt )
  deallocate( sarr,sarr2,ul )
  deallocate( marr,marr2 )
  deallocate( rarr,rarr2 )
  deallocate( zarr,zarr2 )
  deallocate( rrarr,rzarr,zrarr,zzarr,pparr )
  deallocate( aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj )
  deallocate( gaa,gtt,gpp,gat,vthe3,q,ns0,ne0 )

  deallocate(fe_Pijk)
  deallocate(fe_weight)
  deallocate(fe_NN)
  deallocate(fe_dN1)
  deallocate(fe_Eij)
  deallocate(fe_bEij)


  CONTAINS
    !------------------------------------
    SUBROUTINE cmp_phy_qty

     aaa(:)=fem%points(3,ind(:))
     the(:)=fem%points(4,ind(:))
     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,fe_ns)

        ! in (a,t) flux coordinates
        ! fcj: flux coordinate jacobian
        ! Jacobian=(  zaj(1)*rqj(1)-zqj(1)*raj(1)  )
        ! differential area dS = J |grad phi| dadt = J/rj
        fcj = mean( (zaj*rqj-zqj*raj)/rj, fe_ns )

     do i=1,fe_ns
        mt=metric_elements(aaa(i),the(i),0.)
        gaa(i)=metric(1)
        gtt(i)=metric(2)
        gpp(i)=metric(3)
        gat(i)=metric(4)
     end do
        cf1 = mean( gaa, fe_ns )
        cf2 = mean( gat, fe_ns )
        cf3 = mean( gtt, fe_ns )
        cf4 = mean( gpp, fe_ns )

     q(:)=-gFaj(:)/gYaj(:)
     qq = mean( q*q, fe_ns )

     do i=1,fe_ns
        ns0(i)=density(aaa(i),0.0_wp)
        ne0(i)=density(aaa(i),0.0_wp)
     enddo

     do i=1,fe_ns
        vthe3(i)=temperature_e(aaa(i),0.0_wp)
     enddo
     if(ifelectron==1 .or. (ifelectron==2.and.ifdphi==1)) then
        alpha = 0.0_wp
     else
        alpha = abs(csc0) * m0ms / mean(vthe3, fe_ns)
     endif


     !     alpha = 0.0_wp
     !     alpha = alpha/10.0


     beta = - 1./ mean(ne0,fe_ns)
     Zs = 1.
     gs = Zs / mean(Bj*Bj, fe_ns) &
        * mean(ns0, fe_ns) &
        * c0cs_i * msm0_i

    END SUBROUTINE cmp_phy_qty

    FUNCTION mean(vals,nn)
    use precision
    implicit none

    real(wp) :: mean

    integer, intent(in) :: nn
    real(wp), dimension(nn), intent(in) :: vals
    
    mean = sum(vals)/real(nn,wp)
    END FUNCTION mean
end subroutine form_tangent

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
  
  call VecCreateMPIWithArray(mpi_comm_self,1,nn,PETSC_DECIDE,yy,v,ierr)
  !
  call PetscViewerASCIIOpen(mpi_comm_self, str,viewer,ierr)
  call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
  call VecView( v, viewer, ierr )
  call PetscViewerDestroy(viewer,ierr)
  call VecDestroy( v, ierr )
  
end subroutine plot_array
