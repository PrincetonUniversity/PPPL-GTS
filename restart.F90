subroutine restart_write
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use esi_interface
  use diagnosis_array

  implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

  character(len=30) cdum
  character buffer1(MAX_PACKED_LENGTH)
  integer isize,ierr

!  write(cdum,'("RESTART/DATA_RESTART.",i5.5)')mype
  write(cdum,'("RESTART/D",i4.4,"/DATA_RESTART.",i3.3)')myrank_partd,myrank_toroidal
  open(222,file=cdum,status='replace',form='unformatted', access='stream')

  buffer1=''
  isize=pack_sprng(buffer1)   ! get SPRNG seed information

! record particle information for future restart run

  write(222) istep,mi,ntracer,ntracer_energy_ion0,misum
  write(222) zion(1:7,1:mi)
  write(222) zion0(6,1:mi)
  write(222) isize,buffer1
  write(222) sn_box,sn0,sn0_g,sw,rdtem0,rtden,pmarki,markeri,cloud_i    ! these are calculated in load.f90 and should be saved for check point
  write(222) phi  ! --> initial <phi> for poisson solving total phi  (ZF_SEPARATION=0 case)
  write(222) ddT_da  ! need for maintaining constant Ti gradient drive
  write(222) phip00, phi00  ! need for restart when using FEM solver with ZF_SEPARATION=1
  write(222) sumfluxpsi     ! need for restart when using int[I(a)dt] as weighting for enforcing conservation
  if(ifelectron/=0) then
     write(222) me,ntracer_e,ntracer_energy_electron0,mesum
     if(ifdphi/=2) then
        write(222) zeon(1:7,1:me)
     elseif(ifdphi>=2) then
        write(222) zeon(1:8,1:me)     
     endif

     write(222) zeon0(6,1:me)

     write(222) UnNorm_zeta_i(1:mi), UnNorm_theta_i(1:mi)
     write(222) UnNorm_zeta_e(1:mi), UnNorm_theta_e(1:mi)
     write(222) wi_init(1:mi), wi_vpara(1:mi), wi_vperp(1:mi), wi_accel(1:mi), wi_fdamp(1:mi)
     write(222) we_init(1:me), we_vpara(1:me), we_vperp(1:me), we_accel(1:me), we_adiab(1:me), we_fdamp(1:me)

     write(222) pmarke,markere,cloud_e    ! these are calculated in load.f90 and should be saved for check point
!     write(222) dphidt,dphidt00  ! these are saved for check point only when using Jacobian iteration for solving dphidt
  endif

  ! =====================================
  ! MGYOO (07/31/2018)
  ! Add information of multiple tracers
  write(222) num_of_etracers, num_of_itracers, nstep_for_tracer_record
  write(222) etracers_labels(1:num_of_etracers), itracers_labels(1:num_of_itracers)
  write(222) elabels(1:me), ilabels(1:mi)
! =====================================

  close(222)
  
end subroutine restart_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine restart_read
  use global_parameters
  use particle_decomp
  use particle_array
  use field_array
  use esi_interface
  use diagnosis_array

  implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

  SPRNG_POINTER stream
  integer m,isize
  character(len=30) cdum
  character buffer1(MAX_PACKED_LENGTH)
  
!  write(cdum,'("RESTART/DATA_RESTART.",i5.5)')mype
  write(cdum,'("RESTART/D",i4.4,"/DATA_RESTART.",i3.3)')myrank_partd,myrank_toroidal
  open(333,file=cdum,status='old',form='unformatted',access='stream')

! read particle information to restart previous run
  read(333) istart,mi,ntracer,ntracer_energy_ion0,misum
  read(333) zion(1:7,1:mi)
  read(333) zion0(6,1:mi)
  read(333) isize,buffer1
  read(333) sn_box,sn0,sn0_g,sw,rdtem0,rtden,pmarki,markeri,cloud_i
  read(333) phi  ! --> initial <phi> for poisson solving total phi (ZF_SEPARATION=0 case)
  read(333) ddT_da  ! need for maintaining constant Ti gradient drive
  read(333) phip00, phi00  ! need for restart when using FEM solver with ZF_SEPARATION=1
  read(333) sumfluxpsi     ! need for restart when using int[I(a)dt] as weighting for enforcing conservation
  
  istart=istart+1         ! starting step number

  if(ifelectron/=0) then
     read(333) me,ntracer_e,ntracer_energy_electron0,mesum
     if(ifdphi/=2) then
        read(333) zeon(1:7,1:me)
     elseif(ifdphi>=2) then
        read(333) zeon(1:8,1:me)
     endif
     read(333) zeon0(6,1:me)

     read(333) UnNorm_zeta_i(1:mi), UnNorm_theta_i(1:mi)
     read(333) UnNorm_zeta_e(1:mi), UnNorm_theta_e(1:mi)
     read(333) wi_init(1:mi), wi_vpara(1:mi), wi_vperp(1:mi), wi_accel(1:mi), wi_fdamp(1:mi)
     read(333) we_init(1:me), we_vpara(1:me), we_vperp(1:me), we_accel(1:me), we_adiab(1:me), we_fdamp(1:me)

     read(333) pmarke,markere,cloud_e
!     read(333) dphidt,dphidt00  ! these are needed only when using Jacobian iteration for solving dphidt
  endif

  ! =====================================
  ! MGYOO (07/31/2018)
  ! Read information of multiple tracers
  read(333) num_of_etracers, num_of_itracers, nstep_for_tracer_record
  deallocate(etracers_labels,itracers_labels)
  allocate(etracers_labels(num_of_etracers), itracers_labels(num_of_itracers))
  read(333) etracers_labels(1:num_of_etracers), itracers_labels(1:num_of_itracers)
  read(333) elabels(1:me), ilabels(1:mi)
! =====================================

  close(333)

  if(mype==0) print*, 'irun=', irun, '     istart=', istart

  stream=unpack_sprng(buffer1)

  zion0(1:5,1:mi)=zion(1:5,1:mi)
  if(ifelectron/=0) zeon0(1:5,1:mi)=zeon(1:5,1:mi)
!  return

! test domain decomposition
  do m=1,mi
     if(zion(3,m)>zetamax+1.0e-10 .or. zion(3,m)<zetamin-1.0e-10)then
        print *, 'PE=',mype, ' m=',m, ' zion=',zion(3,m)
        stop
     endif
  enddo
  if(ifelectron/=0)then
     do m=1,me
        if(zeon(3,m)>zetamax+1.0e-10 .or. zeon(3,m)<zetamin-1.0e-10)then
           print *, 'PE=',mype, ' m=',m, ' zeon=',zeon(3,m)
           stop
        endif
     enddo
  endif

end subroutine restart_read
