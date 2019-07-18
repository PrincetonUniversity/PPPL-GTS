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
     integer :: order
     integer :: iemax,nmesh_pts
!    integer,dimension(:),pointer :: iet1,iet2,iet3
     integer,dimension(:,:),pointer :: iet
     real(wp),dimension(:,:),pointer :: points
     real(wp),dimension(:,:,:),pointer :: Pijk
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
!     subroutine tri_mesh( this )
!       use fem_data_class
!       implicit none
!       type(fem_data) :: this
!     end subroutine tri_mesh
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
!    deallocate(this%iet1,this%iet2,this%iet3)
     deallocate(this%iet)
     if( associated(this%points) ) deallocate(this%points)
     if( associated(this%Pijk) ) deallocate(this%Pijk)

     this%iemax = 0
!    nullify(this%iet1,this%iet2,this%iet3,this%points)
     nullify(this%iet)
  endif

  ierr = 0
end subroutine delete_fem_data_private

!========================================================================
!                            tri_mesh
!========================================================================    
subroutine tri_mesh( this )
  use precision 
  use global_parameters,only: mpsi,mype,mtheta,zetamin,deltaz,deltat,R_0,pi,numberpe
  use field_array,only: psi_grid,qtinv,q_grid
  !  use petsc_array
  use fem_data_class 
  use esi_interface
  implicit none
  type(fem_data) :: this

  real(wp) zeta, pdum, dt, tdum
  integer i,j,sz,ie,iv
  real(wp), dimension(2) :: aaa,the,rj,zj,Bj
  integer ngvert
  real(wp),dimension(:,:),pointer :: gvert !, tmp
  integer i1,i2,i3,i4,i5,i6

  interface
!    subroutine circ_tri_mesh(npoints,points,order,ngvert,iemax,mype)
!      use precision
!      implicit none
!      integer order, ngvert
!      integer iemax ! in or out
!      integer mype,npoints
!      real(wp) points(5,npoints)
!    end subroutine circ_tri_mesh
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
  call circ_tri_mesh(this%nmesh_pts, this%points, this%order, ngvert, this%iemax, mype)
  
  ! allocate and generate mesh
! allocate( this%iet1(this%iemax), this%iet2(this%iemax), this%iet3(this%iemax) )
  if (this%order .eq. 2) then
  allocate( this%iet(6,this%iemax) )
  else
  allocate( this%iet(3,this%iemax) )
  endif
  !
  ! call mesh_gen( this%iemax, this%iet1, this%iet2, this%iet3, mype )
  call mesh_gen( this%order, this%iemax, this%iet )

  ! extract 2nd order node coordinate
  if (this%order .eq. 2) then
     if(mype==min(1,numberpe)) write(*,*) '	order and global nodes=',this%order, this%nmesh_pts, ngvert
     allocate( gvert(5,ngvert) )
!    allocate( tmp(2,ngvert) )
!     do ie=1,this%iemax
!           i1 = this%iet(1,ie)
!           i2 = this%iet(2,ie)
!           i3 = this%iet(3,ie)
!           gvert(1,i1) = this%points(1,i1)
!           gvert(1,i2) = this%points(1,i2)
!           gvert(1,i3) = this%points(1,i3)
!           gvert(2,i1) = this%points(2,i1)
!           gvert(2,i2) = this%points(2,i2)
!           gvert(2,i3) = this%points(2,i3)
!           i4 = this%iet(4,ie)
!           i5 = this%iet(5,ie)
!           i6 = this%iet(6,ie)
!           gvert(1,i4) = .5*(this%points(1,i1)+this%points(1,i2))
!           gvert(1,i5) = .5*(this%points(1,i2)+this%points(1,i3))
!           gvert(1,i6) = .5*(this%points(1,i3)+this%points(1,i1))
!           gvert(2,i4) = .5*(this%points(2,i1)+this%points(2,i2))
!           gvert(2,i5) = .5*(this%points(2,i2)+this%points(2,i3))
!           gvert(2,i6) = .5*(this%points(2,i3)+this%points(2,i1))
!     enddo
!    call hogrid(this%iemax, this%iet, this%nmesh_pts, this%points(1,:), ngvert, gvert(1,:))
!     do iv=1,ngvert
!       if(abs(gvert(1,iv)-tmp(1,iv)) .gt. 1.e-6) write(*,*) "xcoor is wrong at", iv
!        if( mype.eq.0 ) write(111,111) gvert(1,iv), tmp(1,iv), gvert(1,iv) - tmp(1,iv)
!111 format("xcoor =", 3(1x,d12.6))
!     enddo
!    call hogrid(this%iemax, this%iet, this%nmesh_pts, this%points(2,:), ngvert, gvert(2,:))
!     do iv=1,ngvert
!       if(abs(gvert(2,iv)-tmp(2,iv)) .gt. 1.e-6) write(*,*) "ycoor is wrong at", iv
!        if( mype.eq.0 ) write(112,112) gvert(2,iv), tmp(2,iv), gvert(2,iv) - tmp(2,iv)
!112 format("ycoor =", 3(1x,d12.6))
!     enddo
     do iv=1,5
        call hogrid(this%iemax, this%iet, this%nmesh_pts, this%points(iv,:), ngvert, gvert(iv,:))
     enddo
     DEALLOCATE(this%points)
     this%nmesh_pts =  ngvert
     allocate( this%points(5,this%nmesh_pts) )
     do iv=1,5
        this%points(iv,:) = gvert(iv,:)
     enddo
     DEALLOCATE(gvert)
!     DEALLOCATE(tmp)
  endif

  ! file output -- debug -- change two below also!!!
  if( mype .eq. -1 ) then
     i = 8001+mype
     open(unit=i,file='B2_1.ele')
     write(i,*) this%iemax,' 3 0'
     do ie=1,this%iemax
!       write(i,*) ie,this%iet1(ie),this%iet2(ie),this%iet3(ie)
        write(i,*) ie,this%iet(1,ie),this%iet(2,ie),this%iet(3,ie)
     enddo
     close(i)
  endif

end subroutine tri_mesh

!========================================================================
!                         circ_tri_mesh
!========================================================================
subroutine circ_tri_mesh(npoints,points,order,ngvert,iemax,mype)
  use precision
  use global_parameters,only:mpsi,mtheta
  implicit none
  integer order,ngvert
  integer iemax ! in or out
  integer mype,npoints
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

  call mesh_count(npoints,dpoints,nedges,edges,order,ngvert,iemax,mype)
  DEALLOCATE(dpoints)
  DEALLOCATE(edges)

end subroutine circ_tri_mesh

subroutine hogrid(iemax, iet, insize, inarr, outsize, outarr)
  use precision
  use global_parameters,only: mype
  integer iemax, insize, outsize
  real(wp) inarr(insize)
  real(wp) outarr(outsize)
  integer,dimension(6,iemax) :: iet
  integer:: ie, i1,i2,i3,i4,i5,i6
  integer:: count=0
     do ie=1,iemax
           i1 = iet(1,ie)
           i3 = iet(3,ie)
           i5 = iet(5,ie)
           outarr(i1) = inarr(i1)
           outarr(i3) = inarr(i3)
           outarr(i5) = inarr(i5)
           i2 = iet(2,ie)
           i4 = iet(4,ie)
           i6 = iet(6,ie)
           outarr(i2) = .5*(inarr(i1)+inarr(i3))
           outarr(i4) = .5*(inarr(i3)+inarr(i5))
           outarr(i6) = .5*(inarr(i5)+inarr(i1))
!           if( mype.eq.0 .and. count.lt.1) write(113,113) &
!           ie, inarr(i1), inarr(i3), inarr(i5)
!113 format("xyin =", i7,3(1x,d12.6))
!           if( mype.eq.0 .and. count.lt.1) write(114,114) &
!           ie, outarr(i1), outarr(i3), outarr(i5), outarr(i2), outarr(i4), outarr(i6)
!114 format("xyout =", i7, 6(1x,d12.6))
!           if( mype.eq.0 .and. count.lt.1) write(115,115) &
!           ie, i1, i2, i3, i4, i5, i6
!115 format("ie =", i7, 6(1x,i7))
     enddo
     count=count+1
     if( mype.eq.0 .and. count.eq.mype) write(0,*) '	hogrid is called at', count
end subroutine hogrid

subroutine hogridinv(newsize, newarr, oldsize, oldarr)
  use precision
  integer newsize, oldsize
  real(wp) newarr(newsize)
  real(wp) oldarr(oldsize)
     if(oldsize .lt. newsize) stop '	fem: array size is not right'
     newarr(1:newsize) = oldarr(1:newsize)
end subroutine hogridinv
!  integer iemax, ssize, lsize
!  real(wp) sarr(ssize), larr(lsize)
!  integer,dimension(6,iemax) :: iet
!  integer:: ie, i1,i2,i3
!     if(lsize .lt. ssize) stop '	fem: array size is not right'
!     do ie=1,iemax
!           i1 = iet(1,ie)
!           i2 = iet(2,ie)
!           i3 = iet(3,ie)
!           sarr(i1) = larr(i1)
!           sarr(i2) = larr(i2)
!           sarr(i3) = larr(i3)
!     enddo
