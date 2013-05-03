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
     integer :: numbdy
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
  subroutine new_phit_matrix( this, mat, mindex, bc_type )
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
    integer mindex, bc_type
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
          if(bc_type==1) then      ! for Poisson solver
!cj             if(ZF_SEPARATION==1) then
          ! dirichlet boundary condition at amin & neumann boundary condition at amax
                if( i==0 ) this%gk_matid(it) = -1         ! Dirichlet boundary condiction at inner boundary: phi(i=0)=0
!cj                write(0,*) "	finite element: bc_type=1 ZF_SEPARATION=", ZF_SEPARATION
!cj             else ! ZF_SEPARATION==0                
!cj                if(ifelectron==0)  then ! adiabatic electrons
!cj                   ! NOTE: any combination of following BCs produces a spike at the boundary for ad-electron case
!cj                   ! if( i==0 ) this%gk_matid(it) = -1
!cj                   ! if( i==mpsi ) this%gk_matid(it) = -1  ! Dirichlet boundary condiction at outer boundary: phi(i=mpsi)=0
!cj                else
!cj                   !   if( i==0 ) this%gk_matid(it) = -1
!cj                   !   if( i==mpsi ) this%gk_matid(it) = -1  ! Dirichlet boundary condiction at outer boundary: phi(i=mpsi)=0
!cj                endif
!cj                write(0,*) "	finite element: bc_type=0"
!cj             endif
          else if(bc_type==2) then ! for dphi/dt solver
          ! dirichlet boundary condition at amin & amax
             if( i==0 ) this%gk_matid(it) = -1
             if( i==mpsi ) this%gk_matid(it) = -1
          else if(bc_type==0) then
          ! neumann boundary condition at amin & amax
          else
          ! 3: neumann boundary condition at amin & dirichlet boundary condition at amax
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

     !finite element
     ! order = 0 original solver
     ! order = 1 linear finite element
     ! order = 2 2nd order finite element
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
  public :: nullify_solver,destroy_solver,new_gk_solver, new_fem_solver
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
    call new_phit_matrix( m, Amat, mindex, bc_type )
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
    integer ierr,neq,mindex,mring, i, j, ir, it, ib, ib2, id
    PetscScalar a
    PetscViewer viewer
    integer, dimension(:), allocatable :: bdy
    Vec              vec_yy


    mring=2
    mindex=64*mring+1 ! estimate for matrix size
    neq = sum( mtheta )

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
    call form_tangent(Amat,Mmat,fem%iemax,fem%iet1,fem%iet2,&
         fem%iet3,neq,fem%points,alpha,beta,ierr)
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
#ifdef OLDWAY
             if( bc_type .eq. 1 .and. i==0 ) then
                ! use 0-based bdy index to set Dirichlet bc for lhs matrix
                bdy(ib2)=id
                this%mat%bdy(ib2)=id
                ib2=ib2+1  !count bdy index and assign global index

                ! use 1-based bdy index to set Dirichlet bc for rhs vector
                this%zero01(1+id) = 1.
                this%zero10(1+id) = 0.
             endif
             if( bc_type .eq. 2 .and. (i==0 .or. i==mpsi) ) then
                ! use 0-based bdy index to set Dirichlet bc for lhs matrix
                bdy(ib2)=id
                this%mat%bdy(ib2)=id
                ib2=ib2+1  !count bdy index and assign global index

                ! use 1-based bdy index to set Dirichlet bc for rhs vector
                this%zero01(1+id) = 1.
                this%zero10(1+id) = 0.
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

!   Set Amat=1 on bdy vertex for applying Dirichlet bc
!   using zero-based bdy index
    a = 1.
    call MatZeroRows(Amat,ib,bdy,a,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
    if(mype==0) write(*,*) "	finite elements matrix bdy fix"

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

  end subroutine new_fem_solver


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
    real(wp),intent(inout) :: xx_gk(nn)
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq
    real(doubleprec),allocatable,dimension(:) :: xx, bb 
    integer :: count = 0


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
             xx(ivsl) = xx_gk(ivgts)
          elseif(bb_gk(ivgts).ne.bb(ivsl-itmax+1))then
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

    integer :: count = 5


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
             if( mype.eq.0 ) write(112,112) ir, it, ij, ivgts, rj(1)-R_0, zj(1), xx(ij), rhs(ij)
112 format("init =", 4(i5,2x), 4(1x,d12.6))
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
    call petsc_solve( solver, neq, bb, xx, ierr)

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
    type(solver_data) :: solver
    integer :: nn
    integer :: ierr
    ! locals
    integer :: ir,it,ivgts,ivsl,itq,itmax,ie,iv,jv,kv,neq,mt
    real(doubleprec),allocatable,dimension(:) :: xx, bb 
    real(wp) zeta, pdum, dt, tdum, tmp
    real(wp), dimension(2) :: aaa,the,rj,zj,Bj
    real(wp) cf1, cf2, cf3, cf4, cf5
    real(wp) alpha, beta, ns0, ne0, gs, vthe, Zs
    real(wp), external ::temperature_e, density
    integer, external :: metric_elements

    real(wp) M, L, alfa
    integer :: whichone = 5


    if(mype.eq.0) &
    write(0, *) "	finite element test solve", whichone, "of five"

!cj the following is a combianation of solve and tri_mesh
    neq = sum(mtheta)
    allocate(xx(neq),bb(neq))
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

             M=1.
             L=1.
             alfa=0.

             tdum = real(it)*dt + zeta*qtinv(ir)
             tdum=tdum/(2.0_wp*pi)
             tdum=2.0_wp*pi*(tdum-aint(tdum))

             aaa(1)=pdum
             the(1)=tdum
             call esigetrzb(rj,zj,Bj,aaa,the,1)

             if(whichone==1) then
                !u = (x + y)^2,
                !rhs: c1f*d^2/dx^2 + cf2*d^/dxdy + cf3*d^2/dy^2 + cf4*d/dx + cf5*d/dy = 8 + 4(x+y)
                cf1=1.
                cf2=1.
                cf3=1.
                cf4=1.
                cf5=1.
                xx(ivsl) = ( (rj(1)-R_0) + zj(1) ) * ( (rj(1)-R_0) + zj(1) )
                bb(ivsl) = cf1*2. + cf2*2.*2. + cf3*2. + cf4*2.*((rj(1)-R_0) + zj(1)) + cf5*2.*((rj(1)-R_0) + zj(1))
             else if(whichone==2) then
                !u = sin(pi M x) cos( L y), del u = -( (pi M)^2 + L^2 ) u
                xx(ivsl) = sin(pi*M*(rj(1)-R_0)) * cos(L*zj(1)+alfa)
                bb(ivsl) = (-(pi*M)*(pi*M) - L*L) * xx(ivsl) 
             else if(whichone==3 .or. whichone==4) then
                !u = sin(pi M (a-a0)/(a1-a0)) cos( L theta + alfa),
                !rhs: del u = -( (pi M / (a1-a0))^2 + L^2 ) u
                xx(ivsl) = sin(pi*M*(pdum-a0)/(a1-a0))*cos(L*tdum+alfa)
                bb(ivsl) = (-(pi*M/(a1-a0))*(pi*M/(a1-a0)) - L*L) * xx(ivsl) 

                !whichone=3: in rz coordinates
                !matrix: Amatrix = 
                !      d^2 phi / da ^2  da ^2\     d^2 phi / da  dt     da  dt \   d^2 phi / dt ^2  dt ^2\
                !del u=------- |(--)  +(--)  | + 2 ------- |(--)(--) + (--)(--)| + ------- |(--)  +(--)  |
                !      da^2    \ dx     dy   /     da dt   \ dx  dx     dy  dy /   dt^2    \ dx     dy   /

                !whichone=4: in at coordinates
                !matrix: Amatrix = dnx[i]*dnx[j] + dny[i]*dny[j]
                !      d^2 phi   d^2 phi
                !del u=------- + -------
                !      da^2      dt^2
             else
                !u = sin(pi M (a-a0)/(a1-a0)) cos( L theta + alfa),
                !rhs: del u = -( (pi M / (a1-a0))^2 + L^2 ) u
                xx(ivsl) = sin(pi*M*(pdum-a0)/(a1-a0))*cos(L*tdum+alfa)

                vthe=temperature_e(pdum,0.0_wp)
                alpha = abs(csc0) * m0ms / vthe
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
             if( solver%zero01(ivsl) .gt. 0 ) then
                xx(ivsl) = - xx(ivsl)
             endif
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    
    ! solve GK Poisson equation
    call petsc_solve( solver, neq, bb, xx, ierr)

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

             if(whichone==1) then
                tmp = ( (rj(1)-R_0) + zj(1) ) * ( (rj(1)-R_0) + zj(1) )
             else if(whichone==2) then
                tmp = sin(pi*M*(rj(1)-R_0)) * cos(L*zj(1)+alfa)
             else
                tmp = sin(pi*M*(pdum-a0)/(a1-a0))*cos(L*tdum+alfa)
             endif

             ! print the inner vertex for debugging
             if( mype.eq.0 .and. solver%zero10(ivsl) .gt.0 ) &
             write(111,111) ivsl, xx(ivsl), tmp, abs(tmp-xx(ivsl))/abs(tmp+xx(ivsl))
111 format("finite element test solve=", 1(2x,i7), 3(1x,d12.6))
          end if
       enddo
    enddo
    if(ivsl.ne.neq)stop'solve: ivsl.ne.neq'
    call flush(111)
    if(mype.eq.0) write(0, *) "	finite element test solve end"
    
    deallocate(bb,xx) 

    ! debug
    if( whichone .eq. 0 .and. mype==-10 ) then
       call plot_array( neq, xx, 'x.m' )
       call plot_array( neq, bb, 'b.m' )
       stop 'plot'
    endif

  end subroutine testsolve


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
    PetscScalar,allocatable,dimension(:) :: yy, zz
    PetscViewer viewer
    integer :: m,n,p
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

       ! apply dirichlet bc
       ! bVec = bVec * zero10 + xVec * zero01
       yy=xx
       call VecCreateSeqWithArray(MPI_COMM_SELF,1,nn,yy,vec_yy,ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
       call VecScatterBegin(solver%scat,vec_yy, solver%xVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
       if(ierr .ne. 0) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
       call VecScatterEnd(solver%scat,vec_yy, solver%xVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
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
    if(solver%order .ge. 1) then
       call VecDestroy( rhsVec, ierr )
       call VecDestroy( tmpVec01, ierr )
       call VecDestroy( tmpVec10, ierr )
    endif
  end subroutine petsc_solve
#endif


!========================================================================
!         set up gk_matid for boundary condition
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
          if(bc_type==1) then      ! for Poisson solver
          ! Dirichlet boundary condiction at inner boundary: phi(i=0)=0
             if( i==0 ) this%gk_matid(it) = -1
          else if(bc_type==2) then      ! for dphi/dt solver
             if( i==0 ) this%gk_matid(it) = -1
             if( i==mpsi ) this%gk_matid(it) = -1
          else if(bc_type==3) then
          ! neumann bc at a0, dirichlet bc at a1
             if( i==mpsi ) this%gk_matid(it) = -1
          else                      ! for 
          ! bc_type = 0
          ! neumann bc at both a0 and a1
          endif
!          if(mype .eq. 0 ) write(112,112) it, this%gk_matid(it)
!112 format("new_gk_matid ", 2(2x,i7))
       end do
    end do
    if(id .ne. sum(mtheta)) call MPI_ABORT(mpi_comm_world,1,ierr)
    if(this%gk_matid(sum(mtheta) + mpsi + 2) .ne. -2) call MPI_ABORT(mpi_comm_world,1,ierr)
  end subroutine new_gk_matid

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
  !     points[5,npts]: coordinates of points in cartesian and flux coordinate and 1/q
  !     
  !     output:
  !     Amat: stiffness matrix
  !     Mmat: mass matrix
  !     ierr: error code
  !
  use precision,only:wp,doubleprec
  use particle_decomp,only:myrank_partd
  use global_parameters, only:mype, csc0, m0ms, c0cs_i, msm0_i, ifelectron
  use esi_interface, only : metric
  implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
  integer iemax,npts,ierr
  integer,dimension(:) :: iet1,iet2,iet3
  real(wp) :: points(5,npts)
  real(wp) alpha, beta
  Mat  Amat, Mmat
  real(doubleprec)  xlt(2,3),sarr(3,3),sarr2(3*3),ul(3,4),eps ! therm2d uses *8
  real(doubleprec)           marr(3,3),marr2(3*3)             ! mass matrix
  real(doubleprec)  rrarr(3,3),rzarr(3,3),zrarr(3,3),zzarr(3,3),pparr(3,3),rarr(3,3),zarr(3,3)
  real(doubleprec) c1, c2
  real(doubleprec) :: fcj=1.  ! flux coordinate jacobian
  real(doubleprec) :: cf1=1., cf2=1., cf3=1., cf4=1., cf5=1., cf6=1.
  integer ie,i,j,k,kk,nz,ind2(3),ind3(3),ind(3),nel,nst
  real(wp), dimension(3) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
  real(wp), dimension(3) :: gaa,gtt,gpp,gat,vthe,q,ns0,ne0
  real(wp) :: Zs, gs,qq

  real(wp), external ::temperature_e, density
  integer, external :: metric_elements

  ! Lagrange finite elements
  integer iw, is, il, mt
  integer fe_order, fe_node, fe_ns, fe_nw
  real(doubleprec),dimension(:),allocatable :: fe_weight
  real(doubleprec),dimension(:,:),allocatable :: fe_NN
  real(doubleprec),dimension(:,:,:),allocatable :: fe_dN1, fe_Eij
  integer :: whichone = 5


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

  if(mype .eq. 0) write(0,*) '	form_tangent: iemax=', iemax
  do ie=1,iemax
     ind(1)=iet1(ie)
     ind(2)=iet2(ie)
     ind(3)=iet3(ie)
     if(ind(1) .eq. 0) stop 'ind(1) == 0'
    
     if(whichone.le.3) then
        ! in (r,z) cylindrical coordinates
        xlt(1,:) = points(1,ind(:))
        xlt(2,:) = points(2,ind(:))
     else
        ! in (a,t) flux coordinates
        xlt(1,:) = points(3,ind(:))
        xlt(2,:) = points(4,ind(:))
     endif

     if(whichone == 3) then
        ! in (r,z) cylindrical coordinates
        cf1 = mean( raj*raj + rqj*rqj , 3)
        cf2 = mean( raj*zaj + rqj*zqj , 3)
        cf3 = mean( zaj*zaj + zqj*zqj , 3)
     else if(whichone >= 4) then
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
        call lagrange_element_matrix(fe_ns, fe_nw,  &
                                     xlt(1,:), xlt(2,:),  &
                                     fe_weight, fe_NN, fe_dN1, fe_Eij, &
                                     marr,sarr,rrarr,rzarr,zrarr,zzarr,pparr,rarr,zarr)
#endif

        nz = 0
        do i=1,nel
           k = ind(i)
           nz = nz + 1
           ind3(nz) = i
           ind2(nz) = k - 1 ! zero based
        enddo

        kk = 0
        do i=1,nel
           k = ind(i)
           do j=1,nz
             kk = kk + 1
                if(whichone==1) then
             sarr2(kk) = cf1*rrarr(i,ind3(j)) &
                        +cf2*2.*rzarr(i,ind3(j)) &
                        +cf3*zzarr(i,ind3(j)) &
                        +cf4*rarr(i,ind3(j)) &
                        +cf5*zarr(i,ind3(j))
             marr2(kk) = marr(i,ind3(j))
                else if(whichone==2) then
             sarr2(kk) = sarr(i,ind3(j))
             marr2(kk) = marr(i,ind3(j))
                else if(whichone==3) then
             sarr2(kk) = cf1*rrarr(i,ind3(j)) + 2.*cf2*rzarr(i,ind3(j)) + cf3*zzarr(i,ind3(j))
             marr2(kk) = marr(i,ind3(j))
                else if(whichone==4) then
             sarr2(kk) = fcj*sarr(i,ind3(j))
             marr2(kk) = fcj*marr(i,ind3(j))
                else
             sarr2(kk) = alpha * fcj * marr(i,ind3(j)) &
                       + beta * fcj * gs * ( cf1*rrarr(i,ind3(j)) &
                       +                     cf2*rzarr(i,ind3(j)) &
                       +                     cf2*zrarr(i,ind3(j)) &
                       +                     cf3*zzarr(i,ind3(j)) &
                       +                     cf4*pparr(i,ind3(j))/qq &
                                           )
             marr2(kk) = fcj*marr(i,ind3(j))
                endif
           enddo
        enddo
        !
        call MatSetValues(Amat,nz,ind2,nz,ind2,sarr2,ADD_VALUES,ierr)
        call MatSetValues(Mmat,nz,ind2,nz,ind2,marr2,ADD_VALUES,ierr)
  enddo

  deallocate(fe_weight)
  deallocate(fe_NN)
  deallocate(fe_dN1)
  deallocate(fe_Eij)


  CONTAINS
    !------------------------------------
    SUBROUTINE cmp_phy_qty

     aaa(:)=points(3,ind(:))
     the(:)=points(4,ind(:))
     call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,fe_ns)

        ! in (a,t) flux coordinates
        ! fcj: flux coordinate jacobian
        ! Jacobian=rj(1)*(  zaj(1)*rqj(1)-zqj(1)*raj(1)  )
        fcj = mean( (zaj*rqj-zqj*raj)*rj, fe_ns )

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
        vthe(i)=temperature_e(aaa(i),0.0_wp)
     enddo
     if(ifelectron==1) then
        alpha = 0.0_wp
     else
        alpha = abs(csc0) * m0ms / mean(vthe, fe_ns)
     endif
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
