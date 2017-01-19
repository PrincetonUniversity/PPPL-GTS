subroutine open_files
  use global_parameters
  use particle_decomp
  use field_array
  use diagnosis_array

  implicit none
  
  integer i
  character(len=16) chdum,chdum1,chdum2

  allocate(nf(nbin_a),nf_e(nbin_a))
  do i=1,nbin_a
     nf(i)=50+i
     nf_e(i)=250+i
  end do

  my_nf=120+mype+1
  my_nfL=160+mype+1
  nfile=(mpsi+1)/ncolm
  if((nfile*ncolm).lt.(mpsi+1)) nfile=nfile+1 ! # of pe needed to write ZF
  my_ZL(1)=ncolm*(mype-1)
  my_ZL(2)=min(mpsi,(mype*ncolm-1))

!  if(myrank_toroidal < 10)then
!     write(chdum,'("PHI.00",i1)') myrank_toroidal
!  elseif(myrank_toroidal >= 10 .and. myrank_toroidal < 100)then
!     write(chdum,'("PHI.0",i2)')  myrank_toroidal
!  else
!     write(chdum,'("PHI.",i3)')   myrank_toroidal
!  endif
    write(chdum,'("PHI.",i5.5)')myrank_toroidal
    write(chdum1,'("DEN.",i5.5)')myrank_toroidal
    write(chdum2,'("FLOW.",i5.5)')myrank_toroidal
 
  if(irun==0) then
     if(movie_3D_h5==0 .and. myrank_partd==0) then   ! when movie 3D data output is off
        open(unit=777,file=chdum,form='unformatted')
        write(777) msnap,myrank_toroidal,mpsi,mgrid,mzetamax,ntoroidal,zetamin,zetamax,deltaz,vthc
        write(777) igrid,deltat,qtinv,mtheta,itran,vth_grid
        call flush(777) 
     endif
     if(myrank_partd==1) then   ! myrank_partd==1 will write delta-n
        open(unit=888,file=chdum1,form='unformatted')
!        write(888) msnap,myrank_toroidal,mpsi,mzeta,mgrid,mzetamax,ntoroidal,zetamin,zetamax,deltaz,vthc
        write(888) msnap,myrank_toroidal,mpsi,mgrid,mzetamax,ntoroidal,zetamin,zetamax,deltaz,vthc
        write(888) igrid,deltat,qtinv,mtheta,itran,vth_grid
        call flush(888)
     endif
     if(myrank_partd==2) then   ! myrank_partd==2 will write delta-Vi 
        open(unit=808,file=chdum2,form='unformatted')
!        write(808) msnap,myrank_toroidal,mpsi,mzeta,mgrid,mzetamax,ntoroidal,zetamin,zetamax,deltaz,vthc
        write(808) msnap,myrank_toroidal,mpsi,mgrid,mzetamax,ntoroidal,zetamin,zetamax,deltaz,vthc
        write(808) igrid,deltat,qtinv,mtheta,itran,vth_grid
        call flush(808)
     endif

     if(mype==0)then
      ! Only PE 0 writes into units nf(1:nbin_a). See get_fluxes() routine.
        do i=1,nbin_a
           write(chdum,'("fort.",i2.2)')nf(i)
           open(nf(i),file=chdum,status='replace')
           if(ifelectron/=0) then
              write(chdum,'("fort.",i3.3)')nf_e(i)
              open(nf_e(i),file=chdum,status='replace')
           end if
        end do
     endif
     if(mype.gt.0.and.mype.lt.(nfile+1)) then
        write(chdum,'("fort.",i3.3)')my_nf
        open(my_nf,file=chdum,status='replace')
        write(chdum,'("fort.",i3.3)')my_nfL
        open(my_nfL,file=chdum,status='replace')
     end if

     if(mype==0)open(80,file='fort.80',status='replace',RECL=8192)
     if(mype==0)open(83,file='fort.83',status='replace',RECL=8192)
     if(mype==0)open(800,file='fort.800',status='replace',RECL=8192)
     if(mype==0)open(801,file='fort.801',status='replace',RECL=8192)
     if(mype==numberpe-1)then
        open(90,file='fort.90',status='replace',RECL=8192)
        open(93,file='fort.93',status='replace',RECL=8192)
        open(94,file='fort.94',status='replace',RECL=8192)
        open(95,file='fort.95',status='replace',RECL=8192)
        open(951,file='fort.951',status='replace',RECL=8192)
        open(96,file='fort.96',status='replace',RECL=8192)
        open(97,file='fort.97',status='replace',RECL=8192)
        open(98,file='fort.98',status='replace',RECL=8192)
        open(99,file='fort.99',status='replace',RECL=8192)
        open(85,file='fort.950',status='replace',RECL=8192)
!ZHIXIN add 2 file to record particle/energy flux on grids
        open(86,file='fort.86',status='replace',RECL=8192)
        open(87,file='fort.87',status='replace',RECL=8192)

     end if
     if(mype==1)then
        open(201,file='fort.201',status='replace',RECL=8192)
        open(202,file='fort.202',status='replace',RECL=8192)
        open(203,file='fort.203',status='replace',RECL=8192)
        open(204,file='fort.204',status='replace',RECL=8192)
        open(205,file='fort.205',status='replace',RECL=8192)
        open(206,file='fort.206',status='replace',RECL=8192)
        open(207,file='fort.207',status='replace',RECL=8192)
        open(208,file='fort.208',status='replace',RECL=8192)
        if(ifelectron/=0) then
           open(301,file='fort.301',status='replace',RECL=8192)
           open(302,file='fort.302',status='replace',RECL=8192)
           open(303,file='fort.303',status='replace',RECL=8192)
           open(304,file='fort.304',status='replace',RECL=8192)
           open(305,file='fort.305',status='replace',RECL=8192)
           open(306,file='fort.306',status='replace',RECL=8192)
           open(307,file='fort.307',status='replace',RECL=8192)
           open(308,file='fort.308',status='replace',RECL=8192)
        end if
     end if
     if(mype==numberpe-2)open(91,file='fort.91',status='replace',RECL=8192)

     if(ifelectron/=0.and.mype==numberpe-3)open(193,file='fort.193',status='replace',RECL=8192)
     if(ifelectron/=0.and.mype==numberpe-3)open(194,file='fort.194',status='replace',RECL=8192)
     if(ifelectron/=0.and.mype==numberpe-3)open(195,file='fort.195',status='replace',RECL=8192)
     if(ifelectron/=0.and.mype==numberpe-3)open(199,file='fort.199',status='replace',RECL=8192)
     if(mype==numberpe-4)open(92,file='fort.92',status='replace',RECL=8192)
     if(mype==numberpe-3)open(98,file='fort.98',status='replace',RECL=8192)

  else
     
     if(movie_3D_h5==0 .and. myrank_partd==0) open(unit=777,file=chdum,form='unformatted',position='append')
     if(myrank_partd==1) open(unit=888,file=chdum1,form='unformatted',position='append')
     if(myrank_partd==2) open(unit=808,file=chdum2,form='unformatted',position='append')

     if(mype==0)then
      ! Only PE 0 writes into units nf(1:nbin_a). See get_fluxes() routine.
        do i=1,nbin_a
           write(chdum,'("fort.",i2.2)')nf(i)
           open(nf(i),file=chdum,status='old',position='append')
           if(ifelectron/=0) then 
              write(chdum,'("fort.",i3.3)')nf_e(i)
              open(nf_e(i),file=chdum,status='old',position='append')
           end if
        end do
     endif
     if(mype.gt.0.and.mype.lt.(nfile+1)) then
        write(chdum,'("fort.",i3.3)')my_nf
        open(my_nf,file=chdum,status='old',position='append')
        write(chdum,'("fort.",i3.3)')my_nfL
        open(my_nfL,file=chdum,status='old',position='append')
     end if
     if(mype==0)open(80,file='fort.80',status='old',position='append',RECL=8192)
     if(mype==0)open(83,file='fort.83',status='replace',RECL=8192)
     if(mype==0)open(800,file='fort.800',status='old',position='append',RECL=8192)
     if(mype==0)open(801,file='fort.801',status='old',position='append',RECL=8192)
     if(mype==numberpe-1)then
        open(90,file='fort.90',status='old',position='append',RECL=8192)
        open(93,file='fort.93',status='old',position='append',RECL=8192)
        open(94,file='fort.94',status='old',position='append',RECL=8192)
        open(95,file='fort.95',status='old',position='append',RECL=8192)
        open(951,file='fort.951',status='old',position='append',RECL=8192)
        open(96,file='fort.96',status='old',position='append',RECL=8192)
        open(97,file='fort.97',status='old',position='append',RECL=8192)
        open(99,file='fort.99',status='old',position='append',RECL=8192)
        open(85,file='fort.950',status='old',position='append',RECL=8192)
!ZHIXIN add 2 file to record particle/energy flux on grids
        open(86,file='fort.86',status='old',position='append',RECL=8192)
        open(87,file='fort.87',status='old',position='append',RECL=8192)
     end if
     if(mype==1)then
        open(201,file='fort.201',status='old',position='append',RECL=8192)
        open(202,file='fort.202',status='old',position='append',RECL=8192)
        open(203,file='fort.203',status='old',position='append',RECL=8192)
        open(204,file='fort.204',status='old',position='append',RECL=8192)
        open(205,file='fort.205',status='old',position='append',RECL=8192)
        open(206,file='fort.206',status='old',position='append',RECL=8192)
        open(207,file='fort.207',status='old',position='append',RECL=8192)
        open(208,file='fort.208',status='old',position='append',RECL=8192)
        if(ifelectron/=0) then
           open(301,file='fort.301',status='old',position='append',RECL=8192)
           open(302,file='fort.302',status='old',position='append',RECL=8192)
           open(303,file='fort.303',status='old',position='append',RECL=8192)
           open(304,file='fort.304',status='old',position='append',RECL=8192)
           open(305,file='fort.305',status='old',position='append',RECL=8192)
           open(306,file='fort.306',status='old',position='append',RECL=8192)
           open(307,file='fort.307',status='old',position='append',RECL=8192)
           open(308,file='fort.308',status='old',position='append',RECL=8192)
        end if
     end if
     if(mype==numberpe-2)open(91,file='fort.91',status='old',position='append',RECL=8192)
     if(ifelectron/=0.and.mype==numberpe-3)open(193,file='fort.193',status='old',position='append',RECL=8192)
     if(ifelectron/=0.and.mype==numberpe-3)open(194,file='fort.194',status='old',position='append',RECL=8192)
     if(ifelectron/=0.and.mype==numberpe-3)open(195,file='fort.195',status='old',position='append',RECL=8192)
     if(ifelectron/=0.and.mype==numberpe-3)open(199,file='fort.199',status='old',position='append',RECL=8192)
     if(mype==numberpe-4)open(92,file='fort.92',status='old',position='append',RECL=8192)
     if(mype==numberpe-3)open(98,file='fort.98',status='old',position='append',RECL=8192)
  end if






!    if(myrank_partd==0) then   
!       write(chdum,'("AVW_ion.",i5.5)')myrank_toroidal
!       write(chdum1,'("AVW_eon.",i5.5)')myrank_toroidal
!
!       open(unit=677,file=chdum)
!       open(unit=678,file=chdum1)
!    endif


end subroutine open_files

