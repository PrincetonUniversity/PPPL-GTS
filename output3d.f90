Subroutine output3d
  use precision
  use mpi_defs
  use global_parameters
  use particle_decomp
  use field_array
  use diagnosis_array
  use esi_interface
  use hdf5

  integer i,j,k,kp,n,ij,istatus,ncid,dataid2(3),dimid(2),dataid1(2),one_or_else
  real(wp) dataout(mgrid,2),theta_start(mpsi+1,2),radial(mpsi+1)
  real(wp) r,theta,zeta,theta0
  integer :: rank,ierr
  integer(HSIZE_T), dimension(2) :: count,dimsf
  integer(HSSIZE_T), dimension(2) :: offset
  integer(HSSIZE_T), dimension(7) :: dimsfi = (/0,0,0,0,0,0,0/)
  character(len=1) cdum(3)
  character(len=9) vdum
  character(len=11) fdum

  real(doubleprec) tt0,tt1,rtc

  real(wp), dimension(10) :: aaa,the,rj,zj,Bj

! HDF5 declarations
  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace      ! Dataspace identifier in memory
  integer(HID_T) :: plist_id      ! Property list identifier
  integer(HID_T) :: grp_id        ! Group identifier 
  integer(HID_T) :: time_id       ! Attribute identifier 
  integer(HID_T) :: mpsi_id       ! Attribute identifier 
  integer(HID_T) :: npe_id        ! Attribute identifier 
  integer(HID_T) :: mtheta_id     ! Attribute identifier 
  integer(HID_T) :: radial_id     ! Attribute identifier 
  integer(HID_T) :: itran_id      ! Attribute identifier 
  integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier 
  integer(HID_T) :: a1Dspace_id   ! Attribute Dataspace identifier 
  integer :: ierror               ! Return error flag
  integer :: comm,info
  real(wp)    :: curr_time

  comm = toroidal_comm
  info = MPI_INFO_NULL

! Since we have several MPI processes on each plane, only one of these will
! participate in writing out the data for that plane. Let's pick the processes
! with myrank_partd=0 on each plane.

if(myrank_partd==0)then

! Initializes the HDF5 library and the Fortran90 interface
  call h5open_f(ierror)

  one_or_else=istep/movie_3D_h5

! Write out the grid and its cartesian mapping only once at the beginning.
  if(one_or_else==1) then

   ! Setup file access property list with parallel I/O access.
     call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierror)
     call h5pset_fapl_mpio_f(plist_id,comm,info,ierror)

   ! Filename for grid quantities
     fdum='RUNcoord.h5'

   ! Create the file collectively.
     call h5fcreate_f(fdum,H5F_ACC_TRUNC_F,file_id,ierror,access_prp=plist_id)

   ! Terminates access to a property list
     call h5pclose_f(plist_id,ierror)

   ! Get group id for root group '/'
     call h5gopen_f(file_id,'/',grp_id,ierror)

   ! Create scalar data space for scalar attributes. 
     call h5screate_f(H5S_SCALAR_F,aspace_id,ierror)

   ! Create 1D data space for 1D attributes. 
     rank=1
     dimsf(1)=mpsi+1   ! Dimension of 1D attributes
     call h5screate_simple_f(rank,dimsf,a1Dspace_id,ierror)

   ! Create group datasets.
     call h5dcreate_f(grp_id,"number_of_flux_surfaces",H5T_NATIVE_INTEGER,&
                      aspace_id,mpsi_id,ierror,H5P_DEFAULT_F)
     call h5dcreate_f(grp_id,"number_of_planes",H5T_NATIVE_INTEGER,aspace_id,&
                      npe_id,ierror,H5P_DEFAULT_F)
     call h5dcreate_f(grp_id,"radial_grid",H5T_NATIVE_DOUBLE,a1Dspace_id,&
                      radial_id,ierror,H5P_DEFAULT_F)
     call h5dcreate_f(grp_id,"poloidal_grid",H5T_NATIVE_INTEGER,a1Dspace_id,&
                      mtheta_id,ierror,H5P_DEFAULT_F)
     call h5dcreate_f(grp_id,"index_shift",H5T_NATIVE_INTEGER,a1Dspace_id,&
                      itran_id,ierror,H5P_DEFAULT_F)
     
   ! Write datasets to file. Only the master process writes the datasets.
     if(myrank_toroidal==0)then
        do i=0,mpsi
   !        radial(i+1)=(a0+deltar*real(i))/gyroradius
           radial(i+1)=r_grid(i)/gyroradius
        enddo

        dimsfi(1)=0
        call h5dwrite_f(mpsi_id,H5T_NATIVE_INTEGER,(mpsi+1),dimsfi,ierror)
        call h5dwrite_f(npe_id,H5T_NATIVE_INTEGER,ntoroidal,dimsfi,ierror)
        dimsfi(1)=dimsf(1)
        call h5dwrite_f(radial_id,H5T_NATIVE_DOUBLE,radial,dimsfi,ierror)
        call h5dwrite_f(mtheta_id,H5T_NATIVE_INTEGER,(mtheta+1),dimsfi,ierror)
        call h5dwrite_f(itran_id,H5T_NATIVE_INTEGER,itran,dimsfi,ierror)
     endif

   ! Release datasets ids, datasets data spaces, and group id.
     call h5dclose_f(mpsi_id,ierror)
     call h5dclose_f(npe_id,ierror)
     call h5dclose_f(mtheta_id,ierror)
     call h5dclose_f(radial_id,ierror)
     call h5dclose_f(itran_id,ierror)
     call h5sclose_f(aspace_id,ierror)
     call h5sclose_f(a1Dspace_id,ierror)
     call h5gclose_f(grp_id,ierror)

!  ! Mapping of magnetic coordinates onto cartesian grid and write into file
   ! We explicitely add the zeta=0 poloidal plane since it will be used by
   ! the vizualization software to connect the torus. myrank_toroidal=0 takes
   ! care of this, having one extra plane to write into the output file (kp=1). 
     cdum(1)='X'
     cdum(2)='Y'
     cdum(3)='Z'
     if(myrank_toroidal.eq.0)then
       kp=1   !add the plane zeta=0 to myrank_toroidal=0
     else
       kp=0   !no changes for myrank_toroidal > 0
     endif
     do n=1,3
        do k=1,1+kp
           zeta=zetamin+deltaz*real(k-kp)
           do i=0,mpsi
              r=psi_grid(i)
              theta_start(i+1,k)=zeta*qtinv(i)
              do j=0,mtheta(i)
                 ij=igrid(i)+j
                 theta=deltat(i)*real(j)+zeta*qtinv(i)
                 aaa(1)=r
                 the(1)=theta
                 call esigetrzb(rj,zj,Bj,aaa,the,1)
! grid coordinates (x,y,z), use geometry center as the origin 
                 if(n==1)then
                    dataout(ij,k)=cos(zeta)*rj(1)
                 elseif(n==2)then
                    dataout(ij,k)=-sin(zeta)*rj(1)
                 else
                    dataout(ij,k)=zj(1)
                 endif
              enddo
           enddo
        enddo

!      If n=1, write out the theta_start(mpsi+1,mzeta) array
       if(n.eq.1)then
!        Create the data space for the 2D dataset.
         rank=2
         dimsf(1)=mpsi+1     ! First dimension of theta_start
         dimsf(2)=ntoroidal+1  ! = mzetamax+1 (for the extra plane)
         dimsfi(1)=dimsf(1)
         dimsfi(2)=dimsf(2)
         call h5screate_simple_f(rank,dimsf,filespace,ierror)

!        Create the dataset with default properties.
         call h5dcreate_f(file_id,"theta_start",H5T_NATIVE_DOUBLE,&
                          filespace,dset_id,ierror)
         call h5sclose_f(filespace,ierror)
!        Each process defines dataset in memory and writes it to the hyperslab
!        in the file.
         count(1)=mpsi+1     ! First dimension of theta_start
         offset(1)=0
         if(myrank_toroidal.eq.0)then
           count(2)=2  !2nd dimension (extra plane for myrank_toroidal=0)
           offset(2)=0
         else
           count(2)=1    !2nd dimension of dataout
           offset(2)=myrank_toroidal+1  !shifted by +1 because of extra
         endif                                !plane at myrank_toroidal=0
         call h5screate_simple_f(rank,count,memspace,ierror)

!        Select hyperslab in the file.
         call h5dget_space_f(dset_id,filespace,ierror)
         call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,&
                                    count,ierror)

!        Create property list for collective dataset write
         call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierror)
         call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_INDEPENDENT_F,ierror)

!        Write the dataset collectively.
         call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,theta_start,dimsfi,ierror, &
              file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)

!        Close dataspaces, the dataset and property list.
         call h5sclose_f(filespace,ierror)
         call h5sclose_f(memspace,ierror)
         call h5dclose_f(dset_id,ierror)
         call h5pclose_f(plist_id,ierror)
       endif
 
!      Create the data space for the 2D dataset.
       rank=2
       dimsf(1)=mgrid     ! First dimension of dataout
       dimsf(2)=ntoroidal+1  ! = mzetamax+1 (for the extra plane)
       dimsfi(1)=dimsf(1)
       dimsfi(2)=dimsf(2)
       call h5screate_simple_f(rank,dimsf,filespace,ierror)

!      Create the dataset with default properties.
       call h5dcreate_f(file_id,cdum(n),H5T_NATIVE_DOUBLE,filespace,&
                        dset_id,ierror)
       call h5sclose_f(filespace,ierror)

!      Each process defines dataset in memory and writes it to the hyperslab
!      in the file.
       count(1)=mgrid     ! First dimension of dataout
       offset(1)=0
       if(myrank_toroidal==0)then
         count(2)=2  !2nd dimension (extra plane for myrank_toroidal=0)
         offset(2)=0
       else
         count(2)=1    !2nd dimension of dataout
         offset(2)=myrank_toroidal+1  !shifted by +1 because of extra
       endif                                !plane at myrank_toroidal=0
       call h5screate_simple_f(rank,count,memspace,ierror)

!      Select hyperslab in the file.
       call h5dget_space_f(dset_id,filespace,ierror)
       call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,&
                                  count,ierror)

!      Create property list for collective dataset write
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierror)
!!       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierror)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F , ierror)

       tt0=rtc()
!      Write the dataset collectively.
       call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,dataout, dimsfi,ierror, &
              file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
       tt1=rtc()

!      Close dataspaces, the dataset and property list.
       call h5sclose_f(filespace,ierror)
       call h5sclose_f(memspace,ierror)
       call h5dclose_f(dset_id,ierror)
       call h5pclose_f(plist_id,ierror)

     !!  write(*,*)mype,'  time for collective HDF5 write =',(tt1-tt0)

     enddo

!    Close the file.
     call h5fclose_f(file_id,ierror)
  endif   ! controled by one_or_else

! potential data file
  if(myrank_toroidal==0)then
    kp=1   !add the plane zeta=0 to myrank_toroidal=0
  else
    kp=0   !no changes for myrank_toroidal > 0
  endif
  do k=1,1+kp
     do i=0,mpsi
        do j=0,mtheta(i)
           ij=igrid(i)+j
           dataout(ij,k)=phi(k-kp,ij)
        enddo
     enddo
  enddo

! time label
  n=istep/movie_3D_h5     
  if(n<10)then
     write(fdum,'("PHI_000",i1,".h5")')n
  elseif(n<100)then
     write(fdum,'("PHI_00",i2,".h5")')n
  elseif(n<1000)then
     write(fdum,'("PHI_0",i3,".h5")')n
  else
     write(fdum,'("PHI_",i4,".h5")')n
  endif
 
! Setup file access property list with parallel I/O access.
  call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierror)
  call h5pset_fapl_mpio_f(plist_id,comm,info,ierror)

! Create the file collectively.
  call h5fcreate_f(fdum,H5F_ACC_TRUNC_F,file_id,ierror,access_prp=plist_id)

! Terminates access to a property list
  call h5pclose_f(plist_id,ierror)

! We now write the current calculation time as an attribute to the
! root group. We first need to get the root group id and then create the
! attribute data space, data type, etc. Only one processor needs to write
! the attribute.

! Get group id for root group '/'
  call h5gopen_f(file_id,'/',grp_id,ierror)

! Create scalar data space for the attribute. 
  call h5screate_f(H5S_SCALAR_F,aspace_id,ierror)

! Create dataset with default properties..
  call h5dcreate_f(grp_id,"Time",H5T_NATIVE_DOUBLE,aspace_id,&
                   time_id,ierror,H5P_DEFAULT_F)
     
  curr_time=real(istep)*tstep  ! Attribute data

! Write the time dataset data. Here, dimsfi is just a dummy argument (ignored).
  if(myrank_toroidal.eq.0)then
    call h5dwrite_f(time_id,H5T_NATIVE_DOUBLE,curr_time,dimsfi,ierror)
  endif
     
! Release attribute id, attribute data space, and group id.
  call h5dclose_f(time_id,ierror)
  call h5sclose_f(aspace_id,ierror)
  call h5gclose_f(grp_id,ierror)

! Create the data space for the dataset.
  rank=2
  dimsf(1)=mgrid     ! First dimension of dataout
  dimsf(2)=ntoroidal+1  ! = mzetamax+1
  dimsfi(1)=dimsf(1)
  dimsfi(2)=dimsf(2)
  call h5screate_simple_f(rank,dimsf,filespace,ierror)

! Write dataset name into variable
  write(vdum,'("Potential")')

! Create the dataset with default properties.
  call h5dcreate_f(file_id,vdum,H5T_NATIVE_DOUBLE,filespace,&
                   dset_id,ierror)
  call h5sclose_f(filespace,ierror)

! Each process defines dataset in memory and writes it to the hyperslab
! in the file.
  count(1)=mgrid     ! First dimension of dataout
  offset(1)=0
  if(myrank_toroidal==0)then
    count(2)=2  !2nd dimension (extra plane for myrank_toroidal=0)
    offset(2)=0
  else
    count(2)=1    !2nd dimension of dataout (=mzetamax/ntoroidal)
    offset(2)=myrank_toroidal+1  !shifted by +1 because of extra plane myrank_toroidal=0
  endif
  call h5screate_simple_f(rank,count,memspace,ierror)

! Select hyperslab in the file.
  call h5dget_space_f(dset_id,filespace,ierror)
  call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,&
                             count,ierror)

! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierror)
!!!  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierror)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierror)

  tt0=rtc()
! Write the dataset
  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,dataout, dimsfi,ierror, &
         file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
  tt1=rtc()

! Close dataspaces, the dataset and property list.
  call h5sclose_f(filespace,ierror)
  call h5sclose_f(memspace,ierror)
  call h5dclose_f(dset_id,ierror)
  call h5pclose_f(plist_id,ierror)

! Close the file.
  call h5fclose_f(file_id,ierror)

! Close the HDF5 library and the Fortran90 interface
  call h5close_f(ierror)
  !!   write(*,*)mype,'  time for collective HDF5 Potential write =',(tt1-tt0)

endif   ! if(myrank_partd==0)then

! The processes not participating in writing to file wait at the following
! barrier. We use the partd_comm communicator since it does not require a
! fully global synchronization.
call MPI_BARRIER(partd_comm,ierr)

end

