!========================================================================
!                     new_fem_data_private
!========================================================================  
!
! FEM data class methods --- M. Adams 
!
!===============================================================================
!                       FEM data structures
!===============================================================================
module fem_data_type
  use precision,only:wp
  implicit none
  type fem_data
     integer :: iemax,nmesh_pts
     integer,dimension(:),pointer :: iet1,iet2,iet3
     real(wp),dimension(:,:),pointer :: points
  end type fem_data
end module fem_data_type
!========================================================================
!                        module fem_data_class
!========================================================================
module fem_data_class
  use fem_data_type
  implicit none
  !
contains
  !========================================================================
  !                          new_fem_data
  !========================================================================
  subroutine new_fem_data( this )
    implicit none
    type(fem_data) :: this

    interface
       ! 
!cj       subroutine new_fem_data_private( this, ierr )
!cj         use fem_data_type
!cj         implicit none
!cj         type(fem_data) :: this
!cj         integer ierr
!cj       end subroutine new_fem_data_private
    end interface

    integer ierr

    call new_fem_data_private( this, ierr )

  end subroutine new_fem_data
  !========================================================================
  !                        delete_fem_data
  !========================================================================  
  subroutine delete_fem_data( this )
    implicit none
    type(fem_data) :: this
    ! 
    interface
!cj       subroutine delete_fem_data_private( this, ierr )
!cj         use fem_data_type
!cj         implicit none
!cj         type(fem_data) :: this
!cj         integer ierr
!cj       end subroutine delete_fem_data_private
    end interface
    !
    integer ierr
    call delete_fem_data_private( this, ierr )
  end subroutine delete_fem_data

end module fem_data_class

subroutine new_fem_data_private( this, ierr )
  use precision,only:wp
  use fem_data_class 
  implicit none
  type(fem_data) :: this
  integer,intent(out) :: ierr
  real(wp) xi,xj,xk,yi,yj,yk,bi,bj,bk,ci,cj,ck,delta,fact
  integer i,iv,jv,kv
  !
  interface
     subroutine tri_mesh( this )
       use fem_data_class
       implicit none
       type(fem_data) :: this
     end subroutine tri_mesh
  end interface
  !
  ! form: iet1,iet2,iet3,iemax
  ! triagulate
  ierr = 0

  call tri_mesh( this ) 

end subroutine new_fem_data_private

!========================================================================
!                   subroutine delete_fem_data_private
!========================================================================
subroutine delete_fem_data_private( this, ierr )
  use fem_data_class 
  implicit none
  type(fem_data) :: this
  integer,intent(out) :: ierr
  !
  if( this%iemax > 0 ) then
     deallocate(this%iet1,this%iet2,this%iet3)
     if( associated(this%points) ) deallocate(this%points)

     this%iemax = 0
     nullify(this%iet1,this%iet2,this%iet3,this%points)
  endif

  ierr = 0
end subroutine delete_fem_data_private

!========================================================================
!                            tri_mesh
!========================================================================    
subroutine tri_mesh( this )
  use precision 
  use global_parameters,only: mpsi,mype,mtheta,zetamin,deltaz,deltat,R_0,pi
  use field_array,only: psi_grid,qtinv,q_grid
  !  use petsc_array
  use fem_data_class 
  use esi_interface
  implicit none
  type(fem_data) :: this

  real(wp) zeta, pdum, dt, tdum
  integer i,j,sz,ie,iv
  real(wp), dimension(2) :: aaa,the,rj,zj,Bj

  interface
    subroutine circ_tri_mesh(mtheta,npoints,points,iemax,mype)
      use precision
      implicit none
      integer iemax ! in or out
      integer mype,npoints,mtheta(0:)
      real(wp) points(5,npoints)
    end subroutine circ_tri_mesh
  end interface
  !

  this%nmesh_pts = sum(mtheta)
  allocate( this%points(5,this%nmesh_pts) )

  ! generate vertices: this%points -- grid position in cartitian coordinates (x,y)
  zeta = zetamin + deltaz ! zeta for my plane	
  iv = 0
  do i=0,mpsi
     pdum=psi_grid(i)
     dt = deltat(i)
     sz = mtheta(i) - 1
     do j=0,sz 
        iv= iv + 1 
        tdum = real(j)*dt + zeta*qtinv(i) 
        tdum=tdum/(2.0_wp*pi)
        tdum=2.0_wp*pi*(tdum-aint(tdum))

        aaa(1)=pdum
        the(1)=tdum
        call esigetrzb(rj,zj,Bj,aaa,the,1)
        this%points(1,iv)= rj(1)-R_0   ! --> r; origin of (R.Z) coordinates is at magnetic axis 
        this%points(2,iv)= zj(1)       ! --> z
        this%points(3,iv)= pdum        ! --> a; flux coordinates
        this%points(4,iv)= tdum        ! --> theta
        this%points(5,iv)= q_grid(i)   ! --> q on flux surfaces (0:mpai)
     enddo
  enddo
  if(iv .ne. this%nmesh_pts) stop 'iv .ne. nmesh_pts'

  ! count elements
  call circ_tri_mesh(mtheta, this%nmesh_pts, this%points, this%iemax, mype)
  
  ! allocate and generate mesh
  allocate( this%iet1(this%iemax), this%iet2(this%iemax), this%iet3(this%iemax) )
  !
  call mesh_gen( this%iemax, this%iet1, this%iet2, this%iet3, mype )

  ! file output -- debug -- change two below also!!!
  if( mype .eq. -1 ) then
     i = 8001+mype
     open(unit=i,file='B2_1.ele')
     write(i,*) this%iemax,' 3 0'
     do ie=1,this%iemax
        write(i,*) ie,this%iet1(ie),this%iet2(ie),this%iet3(ie)
     enddo
     close(i)
  endif

end subroutine tri_mesh

!========================================================================
!                         circ_tri_mesh
!========================================================================
subroutine circ_tri_mesh(mtheta,npoints,points,iemax,mype)
  use precision
  use global_parameters,only:mpsi
  implicit none
  integer iemax ! in or out
  integer mype,npoints,mtheta(0:)
  real(wp) points(5,npoints)
  logical just_count
  integer iplot
  integer  nedges,i,j,k,sz,ierr,kk
  integer,allocatable,dimension(:) :: edges
  real(wp),allocatable,dimension(:) :: dpoints

  ALLOCATE(dpoints(2*npoints))
  ALLOCATE(edges(2*npoints))
  
  j = 0
  do i=1,npoints
     j = j + 1
     dpoints(j) = points(1,i)
     j = j + 1
     dpoints(j) = points(2,i)
  enddo

  k = 0
  kk = 0
  do i=0,mpsi
     sz = mtheta(i) - 1
     do j=0,sz
        k = k + 1
        kk = kk + 1
        edges(k) = kk
        k = k + 1
        if( j .eq. sz ) then
           edges(k) = kk - sz
        else
           edges(k) = kk + 1
        endif
     enddo
  enddo
  if(k .ne. 2*npoints) stop 'k .ne. 2*npoints'
  !
  nedges = k/2
  ! file output -- debug
  if( mype .eq. -1 ) then
     i = 8002+mype
     open(unit=i,file='B2_1.edge')
     write(i,*) nedges,' 0'
     k=-1
     do j=1,nedges
        k=k+2
        write(i,*) j, edges(k), edges(k+1)
     enddo
     close(i)
  endif
  
  ! output -- debug
  if( mype .eq. -1 ) then 
     i = 8003+mype
     open(unit=i,file='B2_1.node')
     write(i,*) npoints,' 2 0 0 '
     kk = -1
     do j=1,npoints
        kk = kk + 2
        write(i,*)j,dpoints(kk),dpoints(kk+1)
     enddo
     close(i)
  endif

  call mesh_count(npoints,dpoints,nedges,edges,iemax,mype)
  DEALLOCATE(dpoints)
  DEALLOCATE(edges)

end subroutine circ_tri_mesh
