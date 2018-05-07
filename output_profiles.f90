Subroutine output_profiles
  use precision
  use mpi_defs
  use global_parameters
  use field_array
  use spl_profiles
  use profile_array
  implicit none

  include 'netcdf.inc'


! # of grid in R and Z direction of R-Z plane
  integer,parameter :: nr=300,nz=360,np=1001
  integer i,j,k,ip,ierr,istatus,dataid,ncid,dimid_prof,dimid_mpsi,psiprof_id,&
          Ti_id,Te_id,ni_id,ne_id,a0_id,a1_id,psigrid_id,tstep_id,utime_id,&
          r0_id,ul_id
  real(doubleprec), dimension(np) :: xx,yni,yne,yTi,yTe
  real(doubleprec) thedum,pdum
  real(wp), external :: density,temperature,temperature_e,domega_da

  character(len=12) cdum
  character(len=4) ddum(3)


  thedum=0.0_wp
  do i=1,np
     xx(i)=1.0d0/real(np-1,doubleprec)*real(i-1,doubleprec)
     yni(i)=density(xx(i),thedum)
     yne(i)=yni(i)
     yTi(i)=temperature(xx(i),thedum)
     yTe(i)=temperature_e(xx(i),thedum)
  end do

  if(mype==0)then
! write to data file in NetCDF format
! open netcdf data file
        istatus=nf_create('GTS_backg_profiles.nc',nf_clobber,ncid)
! define data array dimension
        istatus=nf_def_dim(ncid,'profiles_radial_grid',np,dimid_prof)
        istatus=nf_def_dim(ncid,'mpsip1',mpsi+1,dimid_mpsi)
! define data array id
        istatus=nf_def_var(ncid,'tstep',nf_double,0,dimid_prof,tstep_id)
        istatus=nf_def_var(ncid,'time_unit_in_sec',nf_double,0,dimid_prof,utime_id)
        istatus=nf_def_var(ncid,'length_unit_in_cm',nf_double,0,dimid_prof,ul_id)
        istatus=nf_def_var(ncid,'R0',nf_double,0,dimid_prof,r0_id)
        istatus=nf_def_var(ncid,'a0',nf_double,0,dimid_prof,a0_id)
        istatus=nf_def_var(ncid,'a1',nf_double,0,dimid_prof,a1_id)
        istatus=nf_def_var(ncid,'psi_grid',nf_double,1,dimid_mpsi,psigrid_id)
        istatus=nf_def_var(ncid,'psi_grid_profiles',nf_double,1,dimid_prof,psiprof_id)
        istatus=nf_def_var(ncid,'Ti_profile',nf_double,1,dimid_prof,Ti_id)
        istatus=nf_def_var(ncid,'Te_profile',nf_double,1,dimid_prof,Te_id)
        istatus=nf_def_var(ncid,'ni_profile',nf_double,1,dimid_prof,ni_id)
        istatus=nf_def_var(ncid,'ne_profile',nf_double,1,dimid_prof,ne_id)
! end of netcdf definition
        istatus=nf_enddef(ncid)
! write data
        istatus=nf_put_var_double(ncid,tstep_id,tstep)
        istatus=nf_put_var_double(ncid,utime_id,utime)
        istatus=nf_put_var_double(ncid,ul_id,ulength)
        istatus=nf_put_var_double(ncid,r0_id,R_0)
        istatus=nf_put_var_double(ncid,a0_id,a0)
        istatus=nf_put_var_double(ncid,a1_id,a1)
        istatus=nf_put_var_double(ncid,psigrid_id,psi_grid)
        istatus=nf_put_var_double(ncid,psiprof_id,xx)
        istatus=nf_put_var_double(ncid,Ti_id,yTi)
        istatus=nf_put_var_double(ncid,Te_id,yTe)
        istatus=nf_put_var_double(ncid,ni_id,yni)
        istatus=nf_put_var_double(ncid,ne_id,yne)
        istatus=nf_close(ncid)        
  endif

end subroutine output_profiles
