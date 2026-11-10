subroutine diagnosis
   use precision
   use mpi_defs
   use global_parameters
   use particle_decomp
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface

   implicit none

   integer,parameter :: mquantity=20
   integer i,j,k,ij,m,ierror,istatus(MPI_STATUS_SIZE)
   real(wp) r,rdum,q,grada,kappa,etracer,ptracer(4),fdum(17),adum(17),&
      bdum(0:mpsi),tracer(6),tmarker,xnormal(mflux),rfac,dum,&
      epotential,deln2flux(1:mpsi),deln2fluxtmp(1:mpsi)
   character(len=11) cdum

   real(wp),dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

   real(wp), external :: dT_da, temperature

   save etracer,ptracer,xnormal

   ! send tracer information to PE=0
   tracer=0.0

   if(ifelectron==0) then
      if(ntracer>0)then
         tracer(1:5)=zion0(1:5,ntracer)
         tracer(6)=zion(6,ntracer)      ! zion(6,*) is a const. of motion
         if(mype/=0)call MPI_SEND(tracer,6,mpi_Rsize,0,1,MPI_COMM_WORLD,ierror)
      endif
      if(mype==0 .and. ntracer==0)call MPI_RECV(tracer,6,mpi_Rsize,&
         MPI_ANY_SOURCE,1,MPI_COMM_WORLD,istatus,ierror)
   else ! ifelectron=1 or 2
      if(ntracer_e>0) then
         tracer(1:5)=zeon0(1:5,ntracer_e)
         tracer(6)=zeon(6,ntracer_e)
         if(mype/=0) call MPI_SEND(tracer,6,mpi_Rsize,0,1,MPI_COMM_WORLD,ierror)
      endif
      if(mype==0.and.ntracer_e==0) call MPI_RECV(tracer,6,mpi_Rsize,&
         MPI_ANY_SOURCE,1,MPI_COMM_WORLD,istatus,ierror)
   endif

   ! initialize output file
   if((istep-istart+1)==ndiag)then
      if(mype==0)then
         if(irun==0)then

            ! first time run
            open(ihistory,file='history.out',status='replace')
            write(ihistory,101)irun,mquantity,mflux,num_mode,mstep/ndiag
            write(ihistory,102)tstep*real(ndiag)
            mstepall=0

            ! E X B history
            open(444,file='sheareb.out',status='replace')
            write(444,101) ipsi

         else
            ! E X B history
            open(444,file='sheareb.out',status='old',position='append')
            open(444,file='history.out',status='old',position='append')

            ! restart, copy old data to backup file to protect previous run data
            !           open(555,file='history.out',status='old')
            ! # of run
            !           read(555,101)j

            !           write(cdum,'("histry",i1,".bak")') j-10*(j/10)
            !           open(ihistory,file=cdum,status='replace')
            !           write(ihistory,101)j
            !           do i=1,3
            !              read(555,101)j
            !              write(ihistory,101)j
            !           enddo
            ! # of time step
            !           read(555,101)j
            !           write(ihistory,101)j
            !           do i=0,(mquantity+mflux+4*num_mode)*j !i=0 for tstep
            !              read(555,102)dum
            !              write(ihistory,102)dum
            !           enddo
            !           close(555)
            !           close(ihistory)

            ! if everything goes OK, copy old data to history.out
            !           open(666,file=cdum,status='old')
            !           open(ihistory,file='history.out',status='replace')
            !           read(666,101)j
            !           irun=j+1
            !           write(ihistory,101)irun
            !           do i=1,3
            !              read(666,101)j
            !              write(ihistory,101)j
            !           enddo
            !           read(666,101)mstepall
            !           write(ihistory,101) mstep/ndiag
            !           do i=0,(mquantity+mflux+4*num_mode)*mstepall
            !              read(666,102)dum
            !              write(ihistory,102)dum
            !           enddo
            !           close(666)
            !           mstepall=mstepall*ndiag
         endif

         ! tracer particle energy and angular momentum
         aaa(1)=tracer(1)
         the(1)=tracer(2)
         call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
         r=tracer(1)

         !        epotential=gyroradius*r*r*(0.5*flow0+flow1*r/(3.0*a)+0.25*flow2*r*r/(a*a))  ! ????
         epotential=0.0_wp      ! meanless     meanless    meanless    meanless

         if(ifelectron==0) then
            etracer=tracer(6)**2*Bj(1)+0.5*(tracer(4)*Bj(1)*qion)**2/aion+epotential
         else
            etracer=tracer(6)**2*Bj(1)+0.5*(tracer(4)*Bj(1)*qelectron)**2/aelectron+epotential !checked
         end if
         ptracer(1)=tracer(1)
         ptracer(2)=tracer(4)*Fj(1)    ! (Iv_//)/B
         ptracer(3)=-gFaj(1)/gYaj(1)
         ptracer(4)=0.0

         ! flux normalization
         do i=1,mflux
            r=a0+(a1-a0)/real(mflux)*(real(i)-0.5)
            rdum=(r-a0)/djr
            k=int(rdum)
            grada=(real(k+1)-rdum)*flux_it(k)+(rdum-real(k))*flux_it(k+1)   ! =<|grad_a|>

            ! divided by local Kappa_T to obtain chi_i=c_i rho_i
            !           xnormal(i)=1.0/(kappa*kappati*gyroradius)  ! =L_T/rho_i
            rfac=dT_da(r,0.0)/temperature(r,0.0)  ! =dlnT/da
            xnormal(i)=1.0/(abs(rfac)*gyroradius) ! =1/(|dlnT/da|*rho)

            !          xnormal(i)=1.0/(kappa/grada*kappati*gyroradius)   !wx (-dlnT/da)^-1  ????

            write(stdout,*) "kappa_T at radial_bin=",i,rfac,abs(rfac)*grada !=grad_T=L_T^-1

            !           !wx:  kappa/grada*kappati=-dlnT/da (not -dlnT/dr)
         enddo
      endif
      call MPI_BCAST(mstepall,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

   endif

   ! global sum of fluxes
   fdum(1)=efield
   fdum(2)=entropye/vthc**2
   fdum(3)=entropyi
   fdum(4)=dflowi/vthc
   fdum(5)=dflowe/vthc
   fdum(6)=pfluxi/vthc
   fdum(7)=efluxe/vthc**2
   fdum(8:12)=eflux/vthc**3
   fdum(13:17)=rmarker
   call MPI_REDUCE(fdum,adum,17,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
   !  call MPI_REDUCE(fluxpsi,bdum,mpsi+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)

   deln2fluxtmp=0.0
   do i=1,mpsi
      do j=1,mtheta(i)
         ij=igrid(i)+j
         deln2fluxtmp(i)=deln2fluxtmp(i)-phi(k,ij)*phi(k,ij)*evector(2,k,ij)
      enddo
      deln2fluxtmp(i)=deln2fluxtmp(i)/(real(mtheta(i)*mzetamax)*vthc**5)
   enddo
   if(myrank_partd==0)call MPI_REDUCE(deln2fluxtmp,deln2flux,mpsi,&
      mpi_Rsize,MPI_SUM,0,toroidal_comm,ierror)

   if(mype==0)then
      fdum=adum
      !     fluxpsi=bdum

      ! record tracer information
      aaa(1)=tracer(1)
      the(1)=tracer(2)
      call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
      write(ihistory,102) rj(1)-R_0
      write(ihistory,102) zj(1)
      write(ihistory,102) tracer(2)
      write(ihistory,102) tracer(3)
      write(ihistory,102) tracer(5) ! tracer weight

      r=tracer(1)
      q=-gFaj(1)/gYaj(1)
      if(ifelectron==0) then
         write(ihistory,102) tracer(4)*Bj(1)*qion/(aion*vthc)
      else
         write(ihistory,102) tracer(4)*Bj(1)*qelectron/(aelectron*vthc) !=v_//e/v_thi
      end if

      ! energy error

      !     epotential=gyroradius*r*r*(0.5*flow0+flow1*r/(3.0*a)+0.25*flow2*r*r/(a*a))  ! ???
      epotential=0.0_wp      ! meanless     meanless    meanless    meanless

      if(ifelectron==0) then
         write(ihistory,102) (tracer(6)*tracer(6)*Bj(1) &
            +0.5*(tracer(4)*Bj(1)*qion)**2/aion+epotential)/etracer-1.0
      else
         write(ihistory,102) (tracer(6)*tracer(6)*Bj(1) &
            +0.5*(tracer(4)*Bj(1)*qelectron)**2/aelectron+epotential)/etracer-1.0
      end if

      dum=(tracer(1)-ptracer(1))/(0.5*(ptracer(3)+q))-(tracer(4)*Fj(1)-ptracer(2))
      ptracer(1)=tracer(1)
      ptracer(2)=tracer(4)*Fj(1)
      ptracer(3)=q
      ptracer(4)=ptracer(4)+dum
      write(ihistory,102) ptracer(4) ! momentum error

      ! normalization
      fdum(13:17)=max(1.0,fdum(13:17))
      fdum(1)=fdum(1)/real(ntoroidal)
      tmarker=sum(fdum(13:17))
      efluxi=sum(fdum(8:12))/tmarker ! drift speed normalized by V_th
      fdum(2:7)=fdum(2:7)/tmarker
      fdum(8:12)=fdum(8:12)*xnormal(1:5)/fdum(13:17)
      !     fluxpsi(0:mpsi)=fluxpsi(0:mpsi)/real(mzetamax*mtheta(0:mpsi))/vth_grid(0:mpsi)**2
      ! average e*dphi/T on a magnetic surface --- strange average ?

      ! write fluxes to history file
      write(ihistory,102) efluxi,phi(0,igrid(ipsi))/vthc**2,ddeni,ddene,&
         eradial,fdum(1:12)

      ! write mode amplitude history
      do i=1,2
         do j=1,num_mode
            do k=1,2
               write(ihistory,103) amp_mode(k,j,i)
            enddo
         enddo
      enddo

      ! radial-time 2D data
      !     write(444,102)(phip00(i)/vthc,i=1,mpsi,2)
      !     write(444,102)(phi00(i)/vthc**2,i=1,mpsi)
      write(444,102)(fluxpsi(i),i=1,mpsi,2)
      write(444,102)(deln2flux(i),i=1,mpsi,2)
      if(istep==mstep)close(444)
   endif

   101 format(i6)
   102 format(e10.4)
   103 format(e10.4,2x,i4,2x,i4,2x,i4)

end subroutine diagnosis
