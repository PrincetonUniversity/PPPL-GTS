module Ptl_locators_partial

   implicit none

   contains

      subroutine electron_locator_partial(m_start,m_end)
         use precision
         use global_parameters
         use particle_array
         use field_array
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0,m1


         integer m,im,ii,j01,j00
         real(wp) rdum,tdum,delr,delt(0:mpsi),delz,wz1,wp1,     &
            pi2_inv,psitmp,thetatmp,zetatmp,psimax,psimin

         m0=1
         m1=me
         if(present(m_start)) m0=m_start
         if(present(m_end))   m1=m_end


         delr=1.0_wp/deltar
         delt=2.0_wp*pi/deltat
         delz=1.0_wp/deltaz
         pi2_inv=0.5_wp/pi
         psimin=s_grid(0)
         psimax=s_grid(mpsi)


         if(acc_newer_H_zeon) then
            !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif

         !$acc parallel loop present(zeon,jteon0,jteon1,wzeon,wpeon,wteon0,wteon1)
         do m=m0,m1
            !     psitmp=zeon(1,m)      ! =a
            psitmp=zeon(7,m)       ! =s
            thetatmp=zeon(2,m)
            zetatmp=zeon(3,m)

            wzeon(m)=(zetatmp-zetamin)*delz

            rdum=psitmp
            rdum=delr*max(0.0_wp,min(s1-s0,rdum-s0))
            ii=max(0,min(mpsi-1,int(rdum)))
            wp1=rdum-real(ii,wp)
            wpeon(m)=wp1

            ! inner flux surface
            im=ii
            tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
            tdum=(tdum-int(tdum))*delt(im)
            j00=max(0,min(mtheta(im)-1,int(tdum)))
            jteon0(m)=igrid(im)+j00  ! pay attention to j00=0 case
            wteon0(m)=tdum-real(j00,wp)

            ! outer flux surface
            im=ii+1
            tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+100.0_wp
            tdum=(tdum-int(tdum))*delt(im)
            j01=max(0,min(mtheta(im)-1,int(tdum)))
            jteon1(m)=igrid(im)+j01  ! pay attention to j01=0 case
            wteon1(m)=tdum-real(j01,wp)
         enddo

      end subroutine electron_locator_partial

      subroutine ion_locator_at_larmor_points_partial(m_start,m_end)

         use precision
         use global_parameters
         use particle_array
         use field_array
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0,m1

         integer m,i,im,j,k,ip,jt,ii,j11,j10,j01,j00,larmor,ij,ipjt

         real(wp) rdum,tdum,delr,delt(0:mpsi),delz,r,wp1,wp0,wt11,wt10,wt01,wt00,&
            tflr,pi2_inv,psitmp,thetatmp,zetatmp,rhoi,psimax,psimin,pdum,thedum


         integer inner_ip,outter_ip, inner_jt, outter_jt
         real(wp) gyro_phase, intpl_pgyro(larmor_num), intpl_tgyro(larmor_num)


         m0=1
         m1=mi
         if(present(m_start)) m0=m_start
         if(present(m_end))   m1=m_end


         delr=1.0_wp/deltar
         delt=2.0_wp*pi/deltat
         delz=1.0_wp/deltaz
         pi2_inv=0.5_wp/pi
         !  psimin=psi_grid(0)
         !  psimax=psi_grid(mpsi)
         psimin=s_grid(0)
         psimax=s_grid(mpsi)


         if(acc_newer_H_zion) then
            !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_H_zion=.false.
            acc_newer_D_zion=.false.
         endif
         !$acc parallel loop present(zion,jtion0,jtion1,wzion,wpion,wtion0,wtion1,zion_s_on_GyroOrbit,zion_gq_on_GyroOrbit) private(intpl_pgyro,intpl_tgyro)
         do m=m0,m1
            !psitmp=zion(1,m)   ! =a
            psitmp=zion(7,m)   ! =s
            thetatmp=zion(2,m)
            zetatmp=zion(3,m)
            rhoi=zion(6,m)

            r=psitmp
            ip=max(0,min(mpsi,int((r-s0)*delr+0.5)))  ! grid size is uniform in s but not a
            !   jt=max(0,min(mtheta(ip),int(thetatmp*pi2_inv*delt(ip)+0.5)))
            thedum=thetatmp-zetatmp*qtinv(ip)    ! corresponding theta on zeta=0 plane
            thedum=thedum*pi2_inv+100.0_wp
            thedum=2.0_wp*pi*(thedum-int(thedum))       ! bring back to range [0,2pi]
            jt=max(0,min(mtheta(ip),int(thedum*pi2_inv*delt(ip)+0.5)))
            ipjt=igrid(ip)+jt     ! note special case of jt=0
            wzion(m)=(zetatmp-zetamin)*delz


            !==========================================
            ! MGYOO (10/24/2018) : to store information of GC position.
            ! GC's lower-left cell-index(jtion) and linear-interpolation-weights (wpion, wtion0, wtion1)
            inner_ip = max(0,min(mpsi-1,int((r-s0)*delr)))
            thedum=zion(2,m)-zetatmp*qtinv(inner_ip)              ! corresponding theta on zeta=0 plane
            thedum=modulo(thedum,2.0_wp*pi)/deltat(inner_ip)  ! bring back to range [0,2pi] and divde by delta_theta
            inner_jt = int(thedum)
            wtion0(0,m) = thedum - int(thedum)
            jtion0(0,m) = igrid(inner_ip) + inner_jt
            wpion(0,m) =  (zion(7,m) - s_grid(inner_ip))/(s_grid(inner_ip+1)-s_grid(inner_ip))

            outter_ip = inner_ip+1
            thedum=zion(2,m)-zetatmp*qtinv(outter_ip)     ! corresponding theta on zeta=0 plane
            thedum=modulo(thedum,2.0_wp*pi)/deltat(outter_ip)  ! bring back to range [0,2pi] and divde by delta_theta
            outter_jt = int(thedum)
            wtion1(0,m) = thedum - int(thedum)
            jtion1(0,m) = igrid(outter_ip) + outter_jt
            !==========================================


            !==========================================
            ! MGYOO (10/25/2018) : to get linear-interpolation of pgyro and tgyro
            wp0 = 1.0_wp - wpion(0,m)
            wp1 = wpion(0,m)

            wt00 = 1.0_wp - wtion0(0,m)
            wt10 = wtion0(0,m)

            wt01 = 1.0_wp - wtion1(0,m)
            wt11 = wtion1(0,m)

            !$acc loop seq
            do larmor=1, larmor_num
               intpl_pgyro(larmor) = wp0*(wt00*pgyro(larmor,jtion0(0,m)) + wt10*pgyro(larmor,jtion0(0,m)+1)) + wp1*(wt01*pgyro(larmor,jtion1(0,m)) + wt11*pgyro(larmor,jtion1(0,m)+1))
               intpl_tgyro(larmor) = wp0*(wt00*tgyro(larmor,jtion0(0,m)) + wt10*tgyro(larmor,jtion0(0,m)+1)) + wp1*(wt01*tgyro(larmor,jtion1(0,m)) + wt11*tgyro(larmor,jtion1(0,m)+1))
            enddo
            !==========================================

            !$acc loop seq
            do larmor=1, larmor_num
               !  MODIFIED for FINITE GYRORADIUS CASE
               pdum=max(psimin,min(psimax,psitmp+rhoi*intpl_pgyro(larmor)))
               zion_s_on_GyroOrbit(larmor,m)=psitmp+rhoi*intpl_pgyro(larmor)
               rdum=delr*max(0.0_wp,min(s1-s0,pdum-s0))
               ii=max(0,min(mpsi-1,int(rdum)))
               !        wp1=rdum-real(ii)
               !        wpion(larmor,m)=wp1
               wp1=(pdum-s_grid(ii))/(s_grid(ii+1)-s_grid(ii))   !wx  finite rho_i case
               wpion(larmor,m)=wp1                               !wx  finite rho_i case

               ! particle position in theta
               tflr=thetatmp+rhoi*intpl_tgyro(larmor)
               zion_gq_on_GyroOrbit(larmor,m)=modulo(tflr,2.0*pi)

               ! inner flux surface
               im=ii
               tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0_wp
               tdum=(tdum-int(tdum))*delt(im)
               j00=max(0,min(mtheta(im)-1,int(tdum)))
               jtion0(larmor,m)=igrid(im)+j00  ! pay attention to j00=0 case
               wtion0(larmor,m)=tdum-real(j00)

               ! outer flux surface
               im=ii+1
               tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0_wp
               tdum=(tdum-int(tdum))*delt(im)
               j01=max(0,min(mtheta(im)-1,int(tdum)))
               jtion1(larmor,m)=igrid(im)+j01  ! pay attention to j01=0 case
               wtion1(larmor,m)=tdum-real(j01)
            enddo
         enddo

      end subroutine ion_locator_at_larmor_points_partial


end module Ptl_locators_partial

module Ptl_Src_Loss

   use Ptl_locators_partial
   use precision
   implicit none

   ! Following numbers are local things estimated and controlled by each rank.
   ! They are not global things for all ranks (MPI_COMM_WORLD)
   integer :: e_lost_num_BDY0=0, e_lost_num_BDY1=0 ! number of loss particles across boundary per each call of treat_electrons_out_of_boundary
   integer :: e_cum_lost_num_BDY0=0, e_cum_lost_num_BDY1=0 ! cumulative number of lost particles across boundary during a period of time.

   integer :: i_lost_num_BDY0=0, i_lost_num_BDY1=0 ! number of loss particles across boundary at each istep and per each rank
   integer :: i_cum_lost_num_BDY0=0, i_cum_lost_num_BDY1=0 ! Sum of number of loss particles across boundary over all ranks and isteps.

   type FootPrint
      integer :: dim_gq=50, dim_gf=50 ! Default dimension of gq and gf
      real(wp) inv_dgq, inv_dgf
      real(wp) a, s
      integer :: num_eSD, num_iSD
      real(wp),allocatable, dimension(:) :: gq_1D, gf_1D
      real(wp),allocatable, dimension(:,:) :: a_2D, gq_2D, gf_2D, R_2D, Z_2D, BB_2D
      real(wp),allocatable, dimension(:,:,:) :: eScaData, iScaData ! ScalarData(density, temperautre) of electrons and ions
   end type

   type(FootPrint) FP_BDY1, cum_FP_BDY1


   type OutParticles
      integer :: mout ! number of out particles
      integer :: mout_max ! maximum number of out particles
      real(wp) :: Recycling_Ratio=0.7
      real(wp), allocatable, dimension(:,:) :: P, P0
      real(wp), allocatable, dimension(:) :: energy, backT, backf0, labels
   end type


   type(OutParticles) outION

   contains

      subroutine Initialize_Ptl_Src_Loss
         use global_parameters
         implicit none

         if(Flag_Recycling==1) then
            outION%mout_max = mimax/20
            allocate(outION%P(ziondim,outION%mout_max),outION%P0(ziondim,outION%mout_max))
            allocate(outION%energy(outION%mout_max),outION%backT(outION%mout_max),outION%backf0(outION%mout_max),outION%labels(outION%mout_max))
         endif

         call Initialize_FootPrint(FP_BDY1,a1)
         call Initialize_FootPrint(cum_FP_BDY1,a1)

      end subroutine Initialize_Ptl_Src_Loss

      subroutine Initialize_FootPrint(FP,target_a)
         use global_parameters
         implicit none

         type(FootPrint), intent(inout) :: FP
         real(wp),intent(in) :: target_a

         integer Nq,Nf ! Nq=FP%dim_gq, Nf=FP%dim_gf
         integer qq,ff

         FP%dim_gq = 200
         FP%dim_gf = 100

         Nq=FP%dim_gq
         Nf=FP%dim_gf

         FP%inv_dgq = real(Nq-1)/(2.0*pi)
         FP%inv_dgf = real(Nf-1)/(2.0*pi)


         FP%num_eSD = 6 !Marker, Particle Flux, Heat Flux
         FP%num_iSD = 6

         allocate(FP%gq_1D(Nq), FP%gf_1D(Nf))
         allocate(FP%a_2D(Nq,Nf), FP%gq_2D(Nq,Nf), FP%gf_2D(Nq,Nf))
         allocate(FP%R_2D(Nq,Nf), FP%Z_2D(Nq,Nf), FP%BB_2D(Nq,Nf))
         allocate(FP%eScaData(Nq,Nf,FP%num_eSD), FP%iScaData(Nq,Nf,FP%num_iSD))

         FP%eScaData=0.0
         FP%iScaData=0.0

         FP%a_2D = target_a

         do qq=1,Nq
            FP%gq_1D(qq) = (2.0*pi)/(Nq-1)*(qq-1)
            FP%gq_2D(qq,:) = FP%gq_1D(qq)
         enddo

         do ff=1,Nf
            FP%gf_1D(ff) = (2.0*pi)/(Nf-1)*(ff-1)
            FP%gf_2D(:,ff) = FP%gf_1D(ff)
         enddo

         call esigetrzb(FP%R_2D, FP%Z_2D, FP%BB_2D, FP%a_2D, FP%gq_2D, Nq*Nf)

      end subroutine

      subroutine Measure_and_Cumulate_Electron_FootPrint_at_BDY1(FP, energy)
         use global_parameters
         use particle_array
         use field_array,only:q_grid
         use mpi_defs
         implicit none

         type(FootPrint), intent(inout) :: FP
         real(wp), intent(in) :: energy(me)

         real(wp) ratio, wgq0,wgq1,wgf0,wgf1
         real(wp) P_gq, P_gf, P_safety_q
         real(wp) target_gq, target_gf
         integer gqid, gfid, qq, ff

         real(wp) zeon_a_arr(me)
         integer lost_mid_arr
         real(wp), allocatable,dimension(:,:) :: out_zeon, out_zeon0
         real(wp), allocatable,dimension(:) :: delta_f, f_tot, out_energy, out_backTe, out_back_f0_e
         integer m,m_out,vid
         real(wp) time0,time1

         allocate(out_zeon(size(zeon,1),me),out_zeon0(size(zeon0,1),me))
         allocate(out_energy(me),out_backTe(me),out_back_f0_e(me))


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         !$acc update self(zeon(:,:me),zeon0(:,:me),backTe(:me),back_f0_e(:me))

         m_out=0
         do m=1,me
            if(zeon(1,m)>a1) then
               m_out=m_out+1
               out_zeon(:,m_out)  = zeon(:,m)
               out_zeon0(:,m_out) = zeon0(:,m)
               out_energy(m_out) = energy(m)
               out_backTe(m_out) = backTe(m)
               out_back_f0_e(m_out) = back_f0_e(m)
            endif
         enddo

         allocate(delta_f(m_out),f_tot(m_out))


         if(Flag_Load_and_WeightDef==0)then
            delta_f(:) = out_zeon(5,1:m_out)*eNrNm
            f_tot(:)   = eNrNm
         else
            delta_f(:) = out_zeon(5,1:m_out)
            f_tot(:)   = out_zeon(8,1:m_out)
         endif


         do m=1,m_out

            P_gq = modulo(out_zeon(2,m),2.0*pi)
            P_gf = modulo(out_zeon(3,m),2.0*pi)

            P_safety_q = q_grid(mpsi)

            gfid  = (P_gf)*FP%inv_dgf +1 ! last +1 is because of Fortran-style indexing
            wgf1  = (P_gf-FP%gf_1D(gfid))*FP%inv_dgf
            wgf0  = 1.0-wgf1


            ! (At the plane of gfid): the nearst left plane in toroidal direction
            target_gf = FP%gf_1D(gfid)
            target_gq = modulo(P_gq - (P_gf-target_gf)/P_safety_q, 2.0*pi)
            gqid  = (target_gq)*FP%inv_dgq +1
            wgq1  = (target_gq-FP%gq_1D(gqid))*FP%inv_dgq
            wgq0  = 1.0-wgq1

            ratio = wgf0*wgq0; ff=gfid; qq=gqid
            vid=1; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot


            ratio = wgf0*wgq1; ff=gfid; qq=gqid+1
            vid=1; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot


            ! (At the plane of gfid+1): the nearst right plane in toroidal direction
            target_gf = FP%gf_1D(gfid+1)
            target_gq = modulo(P_gq - (P_gf-target_gf)/P_safety_q, 2.0*pi)
            gqid  = (target_gq)*FP%inv_dgq +1
            wgq1  = (target_gq-FP%gq_1D(gqid))*FP%inv_dgq
            wgq0  = 1.0-wgq1

            ratio = wgf1*wgq0; ff=gfid+1; qq=gqid
            vid=1; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot

            ratio = wgf1*wgq1; ff=gfid+1; qq=gqid+1
            vid=1; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backTe(m)*msm0)*delta_f(m) !Temperature definition 1
            vid=3; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%eScaData(qq,ff,vid) = FP%eScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot
         enddo



         deallocate(delta_f,f_tot)
         deallocate(out_zeon,out_zeon0,out_energy,out_backTe,out_back_f0_e)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Electron_FootPrint_BDY1 takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Measure_and_Cumulate_Electron_FootPrint_at_BDY1

      subroutine Measure_and_Cumulate_Ion_FootPrint_at_BDY1(FP, energy)
         use global_parameters
         use particle_array
         use field_array,only:q_grid
         use mpi_defs
         implicit none

         type(FootPrint), intent(inout) :: FP
         real(wp), intent(in) :: energy(mi)

         real(wp) ratio, wgq0,wgq1,wgf0,wgf1
         real(wp) P_gq, P_gf, P_safety_q
         real(wp) target_gq, target_gf
         integer gqid, gfid, qq, ff

         real(wp) zion_a_arr(mi)
         integer lost_mid_arr
         real(wp), allocatable,dimension(:,:) :: out_zion, out_zion0
         real(wp), allocatable,dimension(:) :: delta_f, f_tot, out_energy, out_backT, out_back_f0_i
         integer m,m_out,vid
         real(wp) time0,time1

         allocate(out_zion(size(zion,1),mi),out_zion0(size(zion0,1),mi))
         allocate(out_energy(mi),out_backT(mi),out_back_f0_i(mi))


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         !$acc update self(zion(:,:mi),zion0(:,:mi),backT(:mi),back_f0_i(:mi))

         m_out=0
         do m=1,mi
            if(zion(1,m)>a1) then
               m_out=m_out+1
               out_zion(:,m_out)  = zion(:,m)
               out_zion0(:,m_out) = zion0(:,m)
               out_energy(m_out) = energy(m)
               out_backT(m_out) = backT(m)
               out_back_f0_i(m_out) = back_f0_i(m)
            endif
         enddo

         allocate(delta_f(m_out),f_tot(m_out))


         if(Flag_Load_and_WeightDef==0)then
            delta_f(:) = out_zion(5,1:m_out)*iNrNm
            f_tot(:)   = iNrNm
         else
            delta_f(:) = out_zion(5,1:m_out)
            f_tot(:)   = out_zion0(8,1:m_out)
         endif


         do m=1,m_out

            P_gq = modulo(out_zion(2,m),2.0*pi)
            P_gf = modulo(out_zion(3,m),2.0*pi)

            P_safety_q = q_grid(mpsi)

            gfid  = (P_gf)*FP%inv_dgf +1 ! last +1 is because of Fortran-style indexing
            wgf1  = (P_gf-FP%gf_1D(gfid))*FP%inv_dgf
            wgf0  = 1.0-wgf1


            ! (At the plane of gfid): the nearst left plane in toroidal direction
            target_gf = FP%gf_1D(gfid)
            target_gq = modulo(P_gq - (P_gf-target_gf)/P_safety_q, 2.0*pi)
            gqid  = (target_gq)*FP%inv_dgq +1
            wgq1  = (target_gq-FP%gq_1D(gqid))*FP%inv_dgq
            wgq0  = 1.0-wgq1

            ratio = wgf0*wgq0; ff=gfid; qq=gqid
            vid=1; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backT(m))*delta_f(m) !Temperature definition 1
            vid=3; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot


            ratio = wgf0*wgq1; ff=gfid; qq=gqid+1
            vid=1; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backT(m))*delta_f(m) !Temperature definition 1
            vid=3; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot


            ! (At the plane of gfid+1): the nearst right plane in toroidal direction
            target_gf = FP%gf_1D(gfid+1)
            target_gq = modulo(P_gq - (P_gf-target_gf)/P_safety_q, 2.0*pi)
            gqid  = (target_gq)*FP%inv_dgq +1
            wgq1  = (target_gq-FP%gq_1D(gqid))*FP%inv_dgq
            wgq0  = 1.0-wgq1

            ratio = wgf1*wgq0; ff=gfid+1; qq=gqid
            vid=1; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backT(m))*delta_f(m) !Temperature definition 1
            vid=3; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot

            ratio = wgf1*wgq1; ff=gfid+1; qq=gqid+1
            vid=1; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*delta_f(m) ! Density
            vid=2; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m)-out_backT(m))*delta_f(m) !Temperature definition 1
            vid=3; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*(2.0/3.0*out_energy(m))*delta_f(m) !
            vid=4; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio ! Marker Distribution
            vid=5; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m) ! Integrating f_tot, so n_tot
            vid=6; FP%iScaData(qq,ff,vid) = FP%iScaData(qq,ff,vid) + ratio*f_tot(m)*(2.0/3.0*out_energy(m)) ! Integrating f_tot*energy, so Te_tot
         enddo



         deallocate(delta_f,f_tot)
         deallocate(out_zion,out_zion0,out_energy,out_backT,out_back_f0_i)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Ion_FootPrint_BDY1 takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Measure_and_Cumulate_Ion_FootPrint_at_BDY1

      subroutine Clear_FootPrints_of_Electrons_and_Ions(FP)
         implicit none
         type(FootPrint),intent(inout) :: FP

         FP%eScaData=0
         FP%iScaData=0
      end subroutine

      subroutine MPI_Reduce_FootPrints_of_Electrons_and_Ions(FP,reduce_root)
         use particle_decomp
         use mpi_defs
         implicit none

         type(FootPrint),intent(inout) :: FP
         integer, intent(in) :: reduce_root
         real(wp) tmp_eScaData(FP%dim_gq,FP%dim_gf,FP%num_eSD), tmp_iScaData(FP%dim_gq,FP%dim_gf,FP%num_iSD)
         integer length,ierr

         ! ============================
         ! Electrons
         ! ============================
         tmp_eScaData = FP%eScaData
         length=FP%dim_gq*FP%dim_gf*FP%num_eSD
         call MPI_Reduce(tmp_eScaData,FP%eScaData,length,mpi_Rsize,MPI_SUM,reduce_root,MPI_COMM_WORLD,ierr)

         ! Poloidal (gq) periodic condition
         FP%eScaData(1,:,:) = FP%eScaData(1,:,:) + FP%eScaData(FP%dim_gq,:,:)
         FP%eScaData(FP%dim_gq,:,:) = FP%eScaData(1,:,:)
         ! Toroidal (gf) periodic condition
         FP%eScaData(:,1,:) = FP%eScaData(:,1,:) + FP%eScaData(:,FP%dim_gf,:)
         FP%eScaData(:,FP%dim_gf,:) = FP%eScaData(:,1,:)


         ! ============================
         ! Ions
         ! ============================
         tmp_iScaData = FP%iScaData
         length=FP%dim_gq*FP%dim_gf*FP%num_iSD
         call MPI_Reduce(tmp_iScaData,FP%iScaData,length,mpi_Rsize,MPI_SUM,reduce_root,MPI_COMM_WORLD,ierr)

         ! Poloidal (gq) periodic condition
         FP%iScaData(1,:,:) = FP%iScaData(1,:,:) + FP%iScaData(FP%dim_gq,:,:)
         FP%iScaData(FP%dim_gq,:,:) = FP%iScaData(1,:,:)
         ! Toroidal (gf) periodic condition
         FP%iScaData(:,1,:) = FP%iScaData(:,1,:) + FP%iScaData(:,FP%dim_gf,:)
         FP%iScaData(:,FP%dim_gf,:) = FP%iScaData(:,1,:)

      end subroutine MPI_Reduce_FootPrints_of_Electrons_and_Ions

      subroutine Cumulate_FP_to_cumFP(FP, cum_FP)

         implicit none
         type(FootPrint),intent(inout) :: FP, cum_FP

         cum_FP%eScaData = cum_FP%eScaData + FP%eScaData
         cum_FP%iScaData = cum_FP%iScaData + FP%iScaData

      end subroutine Cumulate_FP_to_cumFP


      subroutine Cumulate_e_lost_num
         implicit none
         e_cum_lost_num_BDY0 = e_cum_lost_num_BDY0 + e_lost_num_BDY0
         e_cum_lost_num_BDY1 = e_cum_lost_num_BDY1 + e_lost_num_BDY1
      end subroutine Cumulate_e_lost_num

      subroutine Cumulate_i_lost_num
         implicit none
         i_cum_lost_num_BDY0 = i_cum_lost_num_BDY0 + i_lost_num_BDY0
         i_cum_lost_num_BDY1 = i_cum_lost_num_BDY1 + i_lost_num_BDY1
      end subroutine Cumulate_i_lost_num

      subroutine Clear_e_cum_lost_num
         implicit none
         e_cum_lost_num_BDY0 = 0
         e_cum_lost_num_BDY1 = 0
      end subroutine Clear_e_cum_lost_num

      subroutine Clear_i_cum_lost_num
         implicit none
         i_cum_lost_num_BDY0 = 0
         i_cum_lost_num_BDY1 = 0
      end subroutine Clear_i_cum_lost_num

      subroutine Clear_OutParticles
         implicit none
         outION%mout = 0
         outION%P = 0.0
         outION%P0 = 0.0
         outION%backT = 0.0
         outION%backf0 = 0.0
         outION%energy = 0.0
         outION%labels = 0.0
      end subroutine

      subroutine Generate_Particles_from_Boundaries_with_LoadBalancing(delta_t_seconds)
         use global_parameters
         use mpi_defs
         implicit none

         real(wp), intent(in) :: delta_t_seconds
         integer ierror
         integer G_e_cum_lost_num_BDY0, G_e_cum_lost_num_BDY1 ! 'G_' stands for the "global" variables
         integer G_i_cum_lost_num_BDY0, G_i_cum_lost_num_BDY1
         integer G_outIon_mout, G_dN_recycling
         integer :: num_of_new_ptls_per_rank=0
         integer :: dNe_BDY0, dNe_BDY1, dNi_BDY0, dNi_BDY1, dN_recycling ! Local (per mpi rank) number of new particles w.r.t ions, electrons, recycling particles
         real(wp) time0, time1

         call MPI_ALLREDUCE(e_cum_lost_num_BDY0, G_e_cum_lost_num_BDY0, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         call MPI_ALLREDUCE(e_cum_lost_num_BDY1, G_e_cum_lost_num_BDY1, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         call MPI_ALLREDUCE(i_cum_lost_num_BDY0, G_i_cum_lost_num_BDY0, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         call MPI_ALLREDUCE(i_cum_lost_num_BDY1, G_i_cum_lost_num_BDY1, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         call MPI_ALLREDUCE(outION%mout, G_outIon_mout, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         dNe_BDY0 = G_e_cum_lost_num_BDY0/numberpe
         dNe_BDY1 = G_e_cum_lost_num_BDY1/numberpe
         dNi_BDY0 = G_i_cum_lost_num_BDY0/numberpe
         dNi_BDY1 = G_i_cum_lost_num_BDY1/numberpe
         if(mype==0) then
            write(stdout,'(A,I8,2(A,I8,A,I8,A),A,I8)') "istep=",istep,  "; G_cum_lost_num  ION(BDY0+BDY1)=(",G_i_cum_lost_num_BDY0," + ", G_i_cum_lost_num_BDY1,")",&
                                                      "; ELECTRON(BDY0+BDY1)=(",G_e_cum_lost_num_BDY0," + ", G_e_cum_lost_num_BDY1,")", "; G_outIon_mout=",G_outIon_mout
            call flush(stdout)
         endif

         if(Flag_Recycling==1) then

#ifdef _MGY_TIMER_
            time0 = MPI_WTIME()
#endif
            dN_recycling =  Produce_Recycling_Ion_and_Electrons_using_OutParticles()
#ifdef _MGY_TIMER_
            time1=MPI_WTIME()
            if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", Produce_Recycling_Ptls  takes", time1-time0, " secs"
            call flush(stdout)
#endif

            call MPI_ALLREDUCE(dN_recycling, G_dN_recycling, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)

            dNi_BDY1 = (G_i_cum_lost_num_BDY1-G_dN_recycling)/numberpe
            dNe_BDY1 = (G_e_cum_lost_num_BDY1-G_dN_recycling)/numberpe

            write(stdout,'(6(A,I8))') "istep=",istep, "; mype=", mype, "; dNi_BDY1 =", dNi_BDY1,"; dNe_BDY1 =", dNe_BDY1 ,"; dN_recycling=", dN_recycling
            call flush(stdout)
         endif



#ifdef _MGY_TIMER_
         time0 = MPI_WTIME()
#endif
         ! Create Src from inner boundary (BDY0)
         call load_electrons_In_or_Out_Flux(dNe_BDY0, delta_t_seconds,Ptl_BC0_val,a0)
         ! Create Src from inner boundary (BDY1)
         call load_electrons_In_or_Out_Flux(dNe_BDY1, delta_t_seconds,Ptl_BC1_val,a1)
         ! Clear accumulated electron lost number
         call Clear_e_cum_lost_num
#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", Inject_electrons_from_Boundaries takes", time1-time0, " secs"
         call flush(stdout)
#endif


#ifdef _MGY_TIMER_
         time0 = MPI_WTIME()
#endif
         ! Create Src from inner boundary (BDY0)
         call load_ions_In_or_Out_Flux(dNi_BDY0, delta_t_seconds,Ptl_BC0_val,a0)
         ! Create Src from inner boundary (BDY1)
         call load_ions_In_or_Out_Flux(dNi_BDY1, delta_t_seconds,Ptl_BC1_val,a1)
         ! Clear accumulated ion lost number
         call Clear_i_cum_lost_num
#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", Inject_ions_from_Boundaries takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Generate_Particles_from_Boundaries_with_LoadBalancing


      subroutine Inject_electrons_from_Boundaries(delta_t_seconds)
         use global_parameters
         use mpi_defs
         implicit none

         real(wp), intent(in) :: delta_t_seconds

         integer ierror
         integer global_e_cum_lost_num_BDY0
         integer global_e_cum_lost_num_BDY1
         integer num_of_new_ptls_per_rank
         real(wp) time0, time1

#ifdef _MGY_TIMER_
         time0 = MPI_WTIME()
#endif
         call MPI_ALLREDUCE(e_cum_lost_num_BDY0, global_e_cum_lost_num_BDY0, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         call MPI_ALLREDUCE(e_cum_lost_num_BDY1, global_e_cum_lost_num_BDY1, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)

         if(mype==0) then
            write(stdout,'(2(A,I8),A,I8,A)') "istep=",istep, "; Inject electrons, global_e_cum_lost_num (BDY0+BDY1)=(",global_e_cum_lost_num_BDY0," + ", global_e_cum_lost_num_BDY1,")"
            call flush(stdout)
         endif

         ! Create Src from inner boundary (BDY1)
         num_of_new_ptls_per_rank = max(0, global_e_cum_lost_num_BDY1/numberpe)
         call load_electrons_In_or_Out_Flux(num_of_new_ptls_per_rank, delta_t_seconds,Ptl_BC1_val,a1)

         ! Create Src from inner boundary (BDY0)
         num_of_new_ptls_per_rank = max(0, global_e_cum_lost_num_BDY0/numberpe)
         call load_electrons_In_or_Out_Flux(num_of_new_ptls_per_rank, delta_t_seconds,Ptl_BC0_val,a0)

         call Clear_e_cum_lost_num

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", Inject_electrons_from_Boundaries takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Inject_electrons_from_Boundaries


      subroutine Inject_ions_from_Boundaries(delta_t_seconds)
         use global_parameters
         use mpi_defs

         implicit none

         real(wp), intent(in) :: delta_t_seconds
         integer ierror
         integer global_i_cum_lost_num_BDY0
         integer global_i_cum_lost_num_BDY1
         integer num_of_new_ptls_per_rank
         real(wp) time0, time1

#ifdef _MGY_TIMER_
         time0 = MPI_WTIME()
#endif
         call MPI_ALLREDUCE(i_cum_lost_num_BDY0, global_i_cum_lost_num_BDY0, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
         call MPI_ALLREDUCE(i_cum_lost_num_BDY1, global_i_cum_lost_num_BDY1, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)

         if(mype==0) then
            write(stdout,'(2(A,I8),A,I8,A)') "istep=",istep, "; Inject ions, global_i_cum_lost_num (BDY0+BDY1)=(",global_i_cum_lost_num_BDY0," + ", global_i_cum_lost_num_BDY1,")"
            call flush(stdout)
         endif

         ! Create Src from inner boundary (BDY1)
         num_of_new_ptls_per_rank = max(0, global_i_cum_lost_num_BDY1/numberpe)
         call load_ions_In_or_Out_Flux(num_of_new_ptls_per_rank, delta_t_seconds,Ptl_BC1_val,a1)

         ! Create Src from inner boundary (BDY0)
         num_of_new_ptls_per_rank = max(0, global_i_cum_lost_num_BDY0/numberpe)
         call load_ions_In_or_Out_Flux(num_of_new_ptls_per_rank, delta_t_seconds,Ptl_BC0_val,a0)

         call Clear_i_cum_lost_num

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", Inject_ions_from_Boundaries takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Inject_ions_from_Boundaries



      subroutine load_electrons_In_or_Out_Flux(num_of_new_ptls, delta_t_seconds, FluxFactor_over_f0, target_a)

         use global_parameters
         use particle_array
         use precision
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         integer,intent(in) :: num_of_new_ptls
         real(wp),intent(in) :: delta_t_seconds
         real(wp),intent(in) :: FluxFactor_over_f0
         real(wp),intent(in) :: target_a

         integer :: Nm,m0,m,tmp_num
         real(wp) inv_Nm
         real(wp) integral_f0
         real(wp) v_abs, v_theta, vpara, vperp, mean_vabs ! v_abs is amplitude of total velocity vector, and v_theta is the angle between vpara and vperp
         real(wp) vmax,vc,c1,c2,c3,c4,f2val,c5,c6,ratio, Rand !For Flag_Load_and_WeightDef==2 case.


         integer :: minpts, ifail=1, ierror
         real(wp) :: acc, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian

         integer sid
         real(wp) P_s,ff,P_q,inv_ds
         real(wp), allocatable ::Te_norm(:), ge_val(:)

         real(wp), allocatable :: wpi(:,:),ga(:),wga(:,:),gat(:),dgr(:),da(:),dgq(:),dgf(:),Er_particle(:)
         !$acc declare create(dgr,da,dgq,dgf,wpi,Er_particle,ga,wga,gat)
         real(wp) fip,r,dela
         real(wp) e1,e2,e3,wz0,wz1,wp0,wt00,tw10,wp1,wt01,wt11,wt10
         real(wp) intp_phi, intp_FLA_Te0, intp_phi00
         real(wp), allocatable, dimension(:) :: AdiabTerm ! q*phi/T , q is negative charge here

         integer ip,ij
         integer gcmotion_eon, gcmotion_eon_em, gcmotion_eon_em_MGY, gcmotion_eon_em_MGY_vec_form



         if(num_of_new_ptls<=0) then
            return
         endif

         m0 = me+1 ! starting index of newly loaded particles

         if(target_a>a0+1e-3 .and. target_a<a1-1e-3) then
            me = me+ num_of_new_ptls
         else
            ! From these temporay particlse, only the num_of_new_ptls particles in the simulation domain will be picked
            me = me+3*num_of_new_ptls
         endif

         allocate(ge_val(m0:me))


         do m=m0,me
            zeon(1,m)=target_a
            zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)
            zeon(7,m)=min(s1,max(s0,zeon(7,m)))

            ! Pick random points (theta,phi) with Helically aligned along B field
            ! zeon(3,m) = 2.0*pi*sprng()
            zeon(3,m) = zetamin + sprng()*(zetamax-zetamin)
            zeon(2,m) = 2.0*pi*sprng()

            elabels(m) = -1.0
         enddo
         !$acc update device(zeon(:,m0:me),elabels(m0:me))

         ! call n_T_u_group_e(m0)   ! get Te, note that back_f0_e cannot be calculated here because it requires zion(4) and zion(6)
         call n_T_u_group_e(m0,me)   ! get Te, note that back_f0_e cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(backTe(m0:me))
         call esigetbf_partial(zeon,m0,me,zeondim)

         allocate(Te_norm(m0:me))
         do m=m0,me
            Te_norm(m) = backTe(m)*msm0
         enddo


         ! ! Do Sampling for the following ranges
         ! ! -umax*sqrt(T/m) < u_para < umax*sqrt(T/m)
         ! ! -umax*sqrt(m*T)/(eB) < rho_para < umax*sqrt(m*T)/(eB)
         ! ! 0 < u_perp < umax*sqrt(2T/m)
         ! ! 0 < sqrt(mu) < sqrt(T/B)*u
         ! zeon(4,m) = (2.0*rand_5d(4,m)-1.0)*umax *sqrt(msm0*Te_norm(m))/(csc0*BB(m)) !neglected electron rotation
         ! zeon(6,m) = sqrt(Te_norm(m)/BB(m))*rand_5d(5,m)*umax

         if(Flag_Load_and_WeightDef==0) then

         elseif(Flag_Load_and_WeightDef==1) then
            do m=m0,me
               v_abs = umax*sqrt(Te_norm(m)*m0ms)*sqrt(sprng()) ! v_r = amax*sqrt(Rand)
               v_theta = pi*sprng() ! For half-circle, theta should be 0<theta<2pi

               ! NOTE : Direction of vpara is determined by the sign of wall_norm_vec*B_cov_a
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)

               !Following line is for the half-circle sampling in the polar coordinates (vpara, vperp)
               ! Marker Distribution Funcion for full range: ge = Nm / ((2pi) * (2pi) * (0.5*pi*v_abs_max^2)) / (Jacobian(x)*2*pi*vperp)
               ! Note that the length in radial direction will be determined later by the each particle's radial velocity
               ! Therefore, this Marker Distribution Function will be actually 3D in space.
               ! NOTE : Since Nm is not determined yet, it will be mulitplied later
               ge_val(m) = 1.0/(4.0*pi*pi*0.5*pi*umax*umax*Te_norm(m)*m0ms)/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)

               zeon(4,m) = msm0*vpara/(csc0*BB(m))
               zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))
            enddo
         elseif(Flag_Load_and_WeightDef==2) then
            do m=m0,me
               ! First, get v_abs using sprng()
               vmax=umax*sqrt(Te_norm(m)*m0ms)
               vc=sqrt(2.0)*sqrt(Te_norm(m)*m0ms)
               c2=0.01
               f2val=0.05
               c1 = -((-2+c2*vmax+f2val*vmax+ c2*vc-f2val*vc)/( vmax*vc))
               c3 = -((2-2*f2val*vmax-c2*vc+ f2val*vc)/( vmax*(vmax-vc)))
               c4 = -((-2+f2val*vmax+c2*vc)/(vmax-vc))
               Rand = sprng()
               if (Rand < c2*vc+c1*vc**2/2) then
                  v_abs = (-c2+sqrt(c2**2+2*c1*Rand))/c1
               else
                  Rand = Rand -(c2*vc+c1*vc**2/2)
                  v_abs = (-c4 + sqrt(2*c3*Rand+(c4+c3*vc)**2))/c3
               endif
               ! Second, get v_theta using sprng()
               ratio=0.15 ! ratio of (max/min)
               c5=1.0/((1.0+ratio)*0.5*pi)
               c6=(1.0-ratio)/((1.0+ratio)*(0.5*pi)**2)
               Rand = sprng()
               v_theta = 0.5*pi + sign(1.0_wp,Rand-0.5)*((c5-sqrt(c5**2-2*c6*abs(Rand-0.5)))/c6);

               ! NOTE : Direction of vpara is determined by the sign of wall_norm_vec*B_cov_a
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)

               ! Marker Distribution Funcion + Normalization with the phase space Jacobians (Jacobian(x)*2*pi*vperp)
               ! Note that the length in radial direction will be determined later by the each particle's radial velocity
               ! Therefore, this Marker Distribution Function will be actually 3D in space.
               if(v_abs<vc) then
                  ge_val(m) = 1.0*(c1*v_abs+c2)*(c5-c6*abs(v_theta-0.5*pi))/(v_abs)/(4.0*pi*pi)/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)
               else
                  ge_val(m) = 1.0*(c3*v_abs+c4)*(c5-c6*abs(v_theta-0.5*pi))/(v_abs)/(4.0*pi*pi)/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)
               endif

               zeon(4,m) = msm0*vpara/(csc0*BB(m))
               zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))
            enddo
         endif
         !$acc update device(zeon(:,m0:me))

         ! call n_T_u_group_e(m0)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))
         call n_T_u_group_e(m0,me,flag_backf0=.true.)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))
         !$acc update self(back_f0_e(m0:me))

         call electron_locator_partial(m0,me) ! This is needed to calculate evector (gradPhi) and its interpolation (wpi)





         allocate(wpi(3,me),ga(me),wga(3,me),gat(me),dgr(me),da(me),dgq(me),dgf(me),Er_particle(me))
         allocate(AdiabTerm(m0:me))



         ! Calculate evector (wpi) first
         !$acc parallel loop present(wzeon,jteon0,wpeon,wteon0,jteon1,wteon1,wpi)
         do m=m0,me
            wz1=wzeon(m)
            wz0=1.0_wp-wz1

            ij=jteon0(m)
            wp0=1.0_wp-wpeon(m)
            wt00=1.0_wp-wteon0(m)
            e1=wp0*wt00*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
            e2=wp0*wt00*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
            e3=wp0*wt00*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            intp_phi = wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0 = wp0*wt00*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=ij+1
            wt10=1.0_wp-wt00
            e1=e1+wp0*wt10*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
            e2=e2+wp0*wt10*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
            e3=e3+wp0*wt10*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            intp_phi=intp_phi+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0=intp_FLA_Te0+wp0*wt10*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=jteon1(m)
            wp1=1.0_wp-wp0
            wt01=1.0_wp-wteon1(m)
            e1=e1+wp1*wt01*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
            e2=e2+wp1*wt01*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
            e3=e3+wp1*wt01*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            intp_phi=intp_phi+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0=intp_FLA_Te0+wp1*wt01*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=ij+1
            wt11=1.0_wp-wt01
            e1=e1+wp1*wt11*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
            e2=e2+wp1*wt11*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
            e3=e3+wp1*wt11*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            intp_phi=intp_phi+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons
            intp_FLA_Te0=intp_FLA_Te0+wp1*wt11*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            wpi(1,m)=e1
            wpi(2,m)=e2
            wpi(3,m)=e3

            r=(zeon(7,m)-s0)/deltar       ! in radial coordinator s
            ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
            fip=r-real(ip,wp)
            intp_phi00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
            if ( Flag_Def_dne_adiab==0 ) then
               AdiabTerm(m)=abs(csc0)*(intp_phi-real(1.0-ZF_SEPARATION,wp)*intp_phi00)/(backTe(m)*msm0) ! Adiabatic term (e*dphi_NZ/Te)
            elseif(Flag_Def_dne_adiab==1) then
               ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
               AdiabTerm(m)=abs(csc0)*(intp_phi+real(ZF_SEPARATION,wp)*intp_phi00)/(backTe(m)*msm0) ! Adiabatic term (e*dphi_tot/Te)
            elseif(Flag_Def_dne_adiab==2) then
               ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) with FLA_Te0
               AdiabTerm(m)=abs(csc0)*(intp_phi+real(ZF_SEPARATION,wp)*intp_phi00)/(intp_FLA_Te0) ! Adiabatic term (e*dphi_tot/Te)
            elseif(Flag_Def_dne_adiab==3) then
               ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) with FLA_Te0
               AdiabTerm(m)=-1.0+exp(abs(csc0)*(intp_phi+real(ZF_SEPARATION,wp)*intp_phi00)/(backTe(m)*msm0)) ! Adiabatic term (e*dphi_tot/Te)
            else
               AdiabTerm(m)=-1e100
               write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
               STOP
            endif
         enddo

         ! Get Er_particle
         Er401_bp=Er401_bp0
         Er_grid=Er_grid0
         dela=1.0_wp/Ydr
         !$acc parallel loop present(zeon,zeon0,Er_particle)
         do m=m0,me
            r=zeon(1,m)-a0
            r=r*dela
            ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
            fip=r-real(ip,wp)
            Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
         enddo

         if(electromagnetic == 1)then
            !$acc parallel loop present(zeon,ga,wga,gat)
            do m=m0,me
               call Get_MP_alpha_amplitude(zeon(1,m),zeon(2,m),zeon(3,m),istep*tstep, ga(m), wga(:,m), gat(m))
            enddo

            if(Flag_GC_Solver==0) then
               m=gcmotion_eon_em_MGY(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,me,zeondim)
            elseif(Flag_GC_Solver==1) then
               m=gcmotion_eon_em_MGY_vec_form(dgr,da,dgq,dgf,zeon,wpi,ga,wga,gat,Er_particle,nonlinear,m0,me,zeondim)
            else
               write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
               call flush(stdout)
               stop
            endif
         else
            m=gcmotion_eon(dgr,da,dgq,dgf,zeon,wpi,Er_particle,nonlinear,me,zeondim)
         endif

         !$acc update self(da(m0:me))

         if(mype==0) then
            write(stdout,'(3(A,ES14.6))') "uveloci=", uveloci,"; utime=", utime, "; multiplied=",utime*uveloci
            write(stdout,'(3(A,ES14.6))') "vpara(m0)/uveloci=", m0ms*csc0*BB(m0)*zeon(4,m0),"; delta_t_seconds=", delta_t_seconds, "; multiplied=",m0ms*csc0*BB(m0)*zeon(4,m0)*uveloci*delta_t_seconds
            write(stdout,'(3(A,ES14.6))') "da(m0)/uveloci=", da(m0),"; delta_t_seconds=", delta_t_seconds, "; multiplied=",da(m0)*uveloci*delta_t_seconds
            call flush(stdout)
         endif


         do m=m0,me
            if(da(m)<0) then
               zeon(1,m)=target_a-1e-5
               zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)
               zeon(7,m)=min(s1,max(s0,zeon(7,m)))
            else
               zeon(1,m)=target_a+1e-5
               zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,acc,dsda,-1,ifail)
               zeon(7,m)=min(s1,max(s0,zeon(7,m)))
            endif
         enddo


         if(target_a<a0+1e-3 .or. target_a>a1-1e-3) then ! If target_a is near the boudnaries, check the particles
            tmp_num=0
            do m=m0,me
               if(zeon(1,m)>a0 .and. zeon(1,m)<a1) then
                  tmp_num = tmp_num+1
               endif
               if(tmp_num>=num_of_new_ptls) then
                  me=m
                  exit
               endif
            enddo
         endif

         tmp_num=me-m0+1
         ! Nm is total number of new particles regarding whole cpu ranks (MPI_COMM_WORLD)
         call MPI_ALLREDUCE(tmp_num,Nm,1,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,ierror)
         inv_Nm=1.0/Nm

         !$acc update device(zeon(:,m0:me),zeon0(:,m0:me))

         !$acc parallel loop present(zeon,zeon0,we_init,we_vpara,we_vperp,we_accel,we_adiab,we_fdamp,UnNorm_theta_e,UnNorm_zeta_e)
         do m=m0,me

            ! Set initial weights of electrons
            ! zeon(5,m) is "df" or "dh"
            ! zeon(8,m) is "f_tot"
            ! zeon(9,m) is "f_AD"
            ! zeon0(9,m) is "f_0"
            ! zeon0(6,m) is "phase volume"

            ! Correction of Marker distribution regarding number of sample particles (Nm)
            ! and the radial length of the sampling space determined by each particle's radial speed
            ge_val(m) = Nm*ge_val(m)/(abs(da(m))*uveloci*1e-2*delta_t_seconds) ! uveloci*1e-2 = unit velocity in (m/s)

            zeon0(6,m) = 1.0/ge_val(m) ! Phase space volume
            integral_f0 = zeon0(6,m)*back_f0_e(m) ! f_0

            ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
            ! zeon0(9,m) is f_0 (@ irk=1)
            zeon(9,m)  = integral_f0*AdiabTerm(m) ! f_AD
            zeon0(9,m) = integral_f0 ! f_0

            ! Note that zeon(8,m)=full_f or full_h
            ! zeon(5,m)=df or dh
            zeon(5,m) = (FluxFactor_over_f0-1.0)*(integral_f0 + int(ifelectron/2)*zeon(9,m)) ! delta_h = (alpha-1)*(f0+f_AD)
            zeon(8,m) = zeon(5,m) + integral_f0 ! full_h = delta_h+f0


            we_init(m)  = zeon(5,m)
            we_vpara(m) = 0
            we_vperp(m) = 0
            we_accel(m) = 0
            we_adiab(m) = 0
            we_fdamp(m) = 0

            zeon0(1:5,m)     = zeon(1:5,m)
            zeon0(7,m)       = zeon(7,m)
            UnNorm_theta_e(m)= zeon(2,m)
            UnNorm_zeta_e(m) = zeon(3,m)
         enddo

         ! Distribute electrons to appropriate MPI ranks
         call treat_electrons_out_of_boundary(m0,me)
         ! call electron_locator_partial(m0,me)
         ! call n_T_u_group_e(m0,me)



         write(stdout,'(6(A,I8))') "istep=",istep, "; mype=",mype,"; m0=",m0,"; me=",me,"; me-m0=",me-m0,"; memax=",memax
         call flush(stdout)
         call MPI_Barrier(MPI_COMM_WORLD,ierror)

         deallocate(Te_norm)
         deallocate(ge_val)
         deallocate(wpi,ga,wga,gat,dgr,da,dgq,dgf,Er_particle)
         deallocate(AdiabTerm)

      end subroutine load_electrons_In_or_Out_Flux

      subroutine load_ions_In_or_Out_Flux(num_of_new_ptls, delta_t_seconds, FluxFactor_over_f0, target_a)

         use global_parameters
         use particle_array
         use precision
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         integer,intent(in) :: num_of_new_ptls
         real(wp),intent(in) :: delta_t_seconds
         real(wp),intent(in) :: FluxFactor_over_f0
         real(wp),intent(in) :: target_a

         integer :: Nm,m0,m,tmp_num
         real(wp) inv_Nm
         real(wp) integral_f0
         real(wp) v_abs, v_theta, vpara, vperp, mean_vabs ! v_abs is amplitude of total velocity vector, and v_theta is the angle between vpara and vperp
         real(wp) vmax,vc,c1,c2,c3,c4,f2val,c5,c6,ratio, Rand !For Flag_Load_and_WeightDef==2 case.


         integer :: minpts, ifail=1, ierror
         real(wp) :: acc, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian

         integer sid
         real(wp) P_s,ff,P_q,inv_ds
         real(wp), allocatable :: gi_val(:)

         real(wp), allocatable :: wpi(:,:),ga(:),wga(:,:),gat(:),dgr(:),da(:),dgq(:),dgf(:),Er_particle(:)
         !$acc declare create(dgr,da,dgq,dgf,wpi,Er_particle,ga,wga,gat)
         real(wp) fip,r,dela
         real(wp) e1,e2,e3,wz0,wz1,wp0,wt00,tw10,wp1,wt01,wt11,wt10
         integer ip,ij,larmor
         integer gcmotion, gcmotion_em, gcmotion_em_MGY, gcmotion_em_MGY_vec_form


         if(num_of_new_ptls<=0) then
            return
         endif

         m0 = mi+1 ! starting index of newly loaded particles

         if(target_a>a0+1e-3 .and. target_a<a1-1e-3) then
            mi = mi+ num_of_new_ptls
         else
            ! From these temporay particlse, only the num_of_new_ptls particles in the simulation domain will be picked
            mi = mi+3*num_of_new_ptls
         endif

         allocate(gi_val(m0:mi))


         do m=m0,mi
            zion(1,m)=target_a
            zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,ifail)
            zion(7,m)=min(s1,max(s0,zion(7,m)))

            ! Pick random points (theta,phi) with Helically aligned along B field
            ! zion(3,m) = 2.0*pi*sprng()
            zion(3,m) = zetamin + sprng()*(zetamax-zetamin)
            zion(2,m) = 2.0*pi*sprng()

            ilabels(m) = -1.0
         enddo
         !$acc update device(zion(:,m0:mi),ilabels(m0:mi))

         ! call n_T_u_group(m0)   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         call n_T_u_group(m0,mi)   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(backu(m0:mi),backT(m0:mi))
         call esigetbf_partial(zion,m0,mi,ziondim)
         backu(m0:mi)=backu(m0:mi)*F(m0:mi)/BB(m0:mi) ! =u_//=I*omega_t/B



         ! ! Do Sampling for the following ranges
         ! ! -umax*sqrt(T/m) < u_para < umax*sqrt(T/m)
         ! ! -umax*sqrt(m*T)/(eB) < rho_para < umax*sqrt(m*T)/(eB)
         ! ! 0 < u_perp < umax*sqrt(2T/m)
         ! ! 0 < sqrt(mu) < sqrt(T/B)*u
         ! zion(4,m) = (2.0*rand_5d(4,m)-1.0)*umax *sqrt(msm0*Te_norm(m))/(csc0*BB(m)) !neglected electron rotation
         ! zion(6,m) = sqrt(backT(m)/BB(m))*rand_5d(5,m)*umax

         if(Flag_Load_and_WeightDef==0) then

         elseif(Flag_Load_and_WeightDef==1) then
            do m=m0,mi
               v_abs = umax*sqrt(backT(m))*sqrt(sprng()) ! v_r = amax*sqrt(Rand)
               v_theta = pi*sprng() ! For half-circle, theta should be 0<theta<2pi

               ! NOTE : Direction of vpara is determined by the sign of wall_norm_vec*B_cov_a
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)

               !Following line is for the half-circle sampling in the polar coordinates (vpara, vperp)
               ! Marker Distribution Funcion for full range: ge = Nm / ((2pi) * (2pi) * (0.5*pi*v_abs_max^2)) / (Jacobian(x)*2*pi*vperp)
               ! Note that the length in radial direction will be determined later by the each particle's radial velocity
               ! Therefore, this Marker Distribution Function will be actually 3D in space.
               ! NOTE : Since Nm is not determined yet, it will be mulitplied later
               gi_val(m) = 1.0/(4.0*pi*pi*0.5*pi*umax*umax*backT(m))/(Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)


               zion(4,m) = vpara/(BB(m))
               zion(4,m) = zion(4,m) + backu(m)/BB(m)
               zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))
            enddo
         elseif(Flag_Load_and_WeightDef==2) then
            do m=m0,mi
               ! First, get v_abs using sprng()
               vmax=umax*sqrt(backT(m))
               vc=sqrt(2.0)*sqrt(backT(m))
               c2=0.01
               f2val=0.05
               c1 = -((-2+c2*vmax+f2val*vmax+ c2*vc-f2val*vc)/( vmax*vc))
               c3 = -((2-2*f2val*vmax-c2*vc+ f2val*vc)/( vmax*(vmax-vc)))
               c4 = -((-2+f2val*vmax+c2*vc)/(vmax-vc))
               Rand = sprng()
               if (Rand < c2*vc+c1*vc**2/2) then
                  v_abs = (-c2+sqrt(c2**2+2*c1*Rand))/c1
               else
                  Rand = Rand -(c2*vc+c1*vc**2/2)
                  v_abs = (-c4 + sqrt(2*c3*Rand+(c4+c3*vc)**2))/c3
               endif
               ! Second, get v_theta using sprng()
               ratio=0.15 ! ratio of (max/min)
               c5=1.0/((1.0+ratio)*0.5*pi)
               c6=(1.0-ratio)/((1.0+ratio)*(0.5*pi)**2)
               Rand = sprng()
               v_theta = 0.5*pi + sign(1.0_wp,Rand-0.5)*((c5-sqrt(c5**2-2*c6*abs(Rand-0.5)))/c6);

               ! NOTE : Direction of vpara is determined by the sign of wall_norm_vec*B_cov_a
               vpara=v_abs*cos(v_theta)
               vperp=v_abs*sin(v_theta)


               ! Marker Distribution Funcion + Normalization with the phase space Jacobians (Jacobian(x)*2*pi*vperp)
               ! Note that the length in radial direction will be determined later by the each particle's radial velocity
               ! Therefore, this Marker Distribution Function will be actually 3D in space.
               if(v_abs<vc) then
                  gi_val(m) = 1.0*(c1*v_abs+c2)*(c5-c6*abs(v_theta-0.5*pi))/(v_abs)/(4.0*pi*pi)/(Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)
               else
                  gi_val(m) = 1.0*(c3*v_abs+c4)*(c5-c6*abs(v_theta-0.5*pi))/(v_abs)/(4.0*pi*pi)/(Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)
               endif


               zion(4,m) = vpara/(BB(m))
               zion(4,m) = zion(4,m) + backu(m)/BB(m)
               zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))
            enddo
         endif
         !$acc update device(zion(:,m0:mi))

         ! call n_T_u_group(m0)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
         call n_T_u_group(m0,mi,flag_backf0=.true.)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
         !$acc update self(back_f0_i(m0:mi))


         call ion_locator_at_larmor_points_partial(m0,mi) ! This is needed to calculate evector (gradPhi) and its interpolation (wpi)



         allocate(wpi(3,mi),ga(mi),wga(3,mi),gat(mi),dgr(mi),da(mi),dgq(mi),dgf(mi),Er_particle(mi))



         ! Calculate evector (wpi) first
         !$acc parallel loop present(zion,zion0)
         do m=m0,mi
            e1=0.0_wp
            e2=0.0_wp
            e3=0.0_wp
            wz1=wzion(m)
            wz0=1.0_wp-wz1

            !$acc loop seq
            do larmor=1, larmor_num
               ij=jtion0(larmor,m)
               wp0=1.0_wp-wpion(larmor,m)
               wt00=1.0_wp-wtion0(larmor,m)
               e1=e1+wp0*wt00*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
               e2=e2+wp0*wt00*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
               e3=e3+wp0*wt00*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

               ij=ij+1
               wt10=1.0_wp-wt00
               e1=e1+wp0*wt10*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
               e2=e2+wp0*wt10*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
               e3=e3+wp0*wt10*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

               ij=jtion1(larmor,m)
               wp1=1.0_wp-wp0
               wt01=1.0_wp-wtion1(larmor,m)
               e1=e1+wp1*wt01*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
               e2=e2+wp1*wt01*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
               e3=e3+wp1*wt01*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))

               ij=ij+1
               wt11=1.0_wp-wt01
               e1=e1+wp1*wt11*(wz0*evector(1,0,ij)+wz1*evector(1,1,ij))
               e2=e2+wp1*wt11*(wz0*evector(2,0,ij)+wz1*evector(2,1,ij))
               e3=e3+wp1*wt11*(wz0*evector(3,0,ij)+wz1*evector(3,1,ij))
            enddo

            wpi(1,m)=e1*inv_larmor_num
            wpi(2,m)=e2*inv_larmor_num
            wpi(3,m)=e3*inv_larmor_num
         enddo

         ! Get Er_particle
         Er401_bp=Er401_bp0
         Er_grid=Er_grid0
         dela=1.0_wp/Ydr
         !$acc parallel loop present(zion,zion0,Er_particle)
         do m=m0,mi
            r=zion(1,m)-a0
            r=r*dela
            ip=max(0,min(je-1,int(r)))      ! ip [0,je-1]
            fip=r-real(ip,wp)
            Er_particle(m)=fip*Er401_bp(ip+1)+(1.0_wp-fip)*Er401_bp(ip) ! Er at particle's position
         enddo

         if(electromagnetic == 1)then
            !$acc parallel loop present(zion,ga,wga,gat)
            do m=m0,mi
               call Get_MP_alpha_amplitude(zion(1,m),zion(2,m),zion(3,m),istep*tstep, ga(m), wga(:,m), gat(m))
            enddo

            if(Flag_GC_Solver==0) then
               m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,wpi,ga,wga,gat,Er_particle,nonlinear,mi,ziondim)
            elseif(Flag_GC_Solver==1) then
               m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,wpi,ga,wga,gat,Er_particle,nonlinear,m0,mi,ziondim)
            else
               write(stdout,*) "Wrong Flag_GC_Solver=", Flag_GC_Solver
               call flush(stdout)
               stop
            endif
         else
            m=gcmotion(dgr,da,dgq,dgf,zion,wpi,Er_particle,nonlinear,mi,ziondim)
         endif

         !$acc update self(da(m0:mi))

         if(mype==0) then
            write(stdout,'(3(A,ES14.6))') "uveloci=", uveloci,"; utime=", utime, "; multiplied=",utime*uveloci
            write(stdout,'(3(A,ES14.6))') "vpara(m0)/uveloci=", BB(m0)*zion(4,m0),"; delta_t_seconds=", delta_t_seconds, "; multiplied=",m0ms*csc0*BB(m0)*zion(4,m0)*uveloci*delta_t_seconds
            write(stdout,'(3(A,ES14.6))') "da(m0)/uveloci=", da(m0),"; delta_t_seconds=", delta_t_seconds, "; multiplied=",da(m0)*uveloci*delta_t_seconds
            call flush(stdout)
         endif

         ! Set initial weights of electrons
         ! zion(5,m) is "df" or "dh"
         ! zion0(8,m) is "f_tot"

         do m=m0,mi
            if(da(m)<0) then
               zion(1,m)=target_a-1e-5
               zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,ifail)
               zion(7,m)=min(s1,max(s0,zion(7,m)))
            else
               zion(1,m)=target_a+1e-5
               zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,acc,dsda,-1,ifail)
               zion(7,m)=min(s1,max(s0,zion(7,m)))
            endif
         enddo


         if(target_a<a0+1e-3 .or. target_a>a1-1e-3) then
            tmp_num=0
            do m=m0,mi
               if(zion(1,m)>a0 .and. zion(1,m)<a1) then
                  tmp_num = tmp_num+1
               endif
               if(tmp_num>=num_of_new_ptls) then
                  mi=m
                  exit
               endif
            enddo
         endif

         tmp_num=mi-m0+1
         ! Nm is total number of new particles regarding whole cpu ranks (MPI_COMM_WORLD)
         call MPI_ALLREDUCE(tmp_num,Nm,1,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD,ierror)
         inv_Nm=1.0/Nm

         !$acc update device (zion(:,m0:mi),zion0(:,m0:mi))

         !$acc parallel loop present(zion,zion0,wi_init,wi_vpara,wi_vperp,wi_accel,wi_fdamp,UnNorm_theta_i,UnNorm_zeta_i)
         do m=m0,mi

            ! Set initial weights of ions
            ! zion(5,m) is "df"
            ! zion0(8,m) is "f_tot"
            ! zion0(6,m) is "phase volume"

            ! Correction of Marker distribution regarding number of sample particles (Nm)
            ! and the radial length of the sampling space determined by each particle's radial speed
            gi_val(m) = Nm*gi_val(m)/(abs(da(m))*uveloci*1e-2*delta_t_seconds) ! uveloci*1e-2 = unit velocity in (m/s)


            zion0(6,m) = 1.0/gi_val(m) ! Phase space volume
            integral_f0 = zion0(6,m)*back_f0_i(m) !f_0

            zion0(8,m) = FluxFactor_over_f0*integral_f0 ! full_f
            zion(5,m)  = zion0(8,m)-integral_f0 ! df


            wi_init(m)  = zion(5,m)
            wi_vpara(m) = 0
            wi_vperp(m) = 0
            wi_accel(m) = 0
            wi_fdamp(m) = 0

            zion0(1:5,m)     = zion(1:5,m)
            zion0(7,m)       = zion(7,m)
            UnNorm_theta_i(m)= zion(2,m)
            UnNorm_zeta_i(m) = zion(3,m)
         enddo




         ! Distribute electrons to appropriate MPI ranks
         call treat_ions_out_of_boundary(m0,mi)
         ! call ion_locator_at_larmor_points_partial(m0,mi)
         ! call n_T_u_group(m0,mi)


         write(stdout,'(6(A,I8))') "istep=",istep, "; mype=",mype,"; m0=",m0,"; mi=",mi,"; mi-m0=",mi-m0,"; mimax=",mimax
         call flush(stdout)
         call MPI_Barrier(MPI_COMM_WORLD,ierror)

         deallocate(gi_val)
         deallocate(wpi,ga,wga,gat,dgr,da,dgq,dgf,Er_particle)

      end subroutine load_ions_In_or_Out_Flux


      subroutine treat_electrons_out_of_boundary(m_start, m_end)

         use global_parameters
         use particle_array
         use field_array
         use function_gg_module
         use Ptl_locators_partial
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0,m1
         integer m,i,j
         real(wp) psimax,psimin
         integer lost_mid_arr(memax) ! index of loss particles
         integer lost_num ! number of loss particles across boundary
         integer, allocatable :: fill_mid_arr(:)
         integer fill_num
         integer mtop, lasth
         integer dim_zeon, dim_zeon0, dim_RK_cum_zeon
         real(wp) zeon_a_arr(me)

         m0=1
         m1=me
         if(present(m_start)) m0=m_start
         if(present(m_end))   m1=m_end


         if(acc_newer_H_zeon) then
            !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif

         dim_zeon = size(zeon,1)
         dim_zeon0 = size(zeon0,1)
         dim_RK_cum_zeon = size(RK_cum_zeon,1)

         psimax=psi_grid(mpsi)
         psimin=psi_grid(0)

         if(Flag_Ptl_BC0_type==0 .or. Flag_Ptl_BC1_type==0) then
            !$acc parallel loop present(zeon,zeon0)
            do m=m0,m1
               if(Flag_Ptl_BC1_type==0 .and. zeon(1,m) > psimax)then
                  zeon(1,m)=psimax
                  zeon(2,m)=2.0_wp*pi-zeon0(2,m)
                  zeon(3,m)=zeon0(3,m)
                  zeon(4,m)=zeon0(4,m)
                  zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
                  zeon(7,m)=zeon0(7,m)

                  if(Flag_w_Evol==3) zeon(8,m)=zeon0(8,m)
               elseif(Flag_Ptl_BC0_type==0 .and. zeon(1,m) < psimin)then
                  zeon(1,m)=psimin
                  zeon(2,m)=2.0_wp*pi-zeon0(2,m)
                  zeon(3,m)=zeon0(3,m)
                  zeon(4,m)=zeon0(4,m)
                  zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
                  zeon(7,m)=zeon0(7,m)
               endif
            enddo
            acc_newer_D_zeon=.true.

         endif

         if(Flag_Ptl_BC0_type==1 .or. Flag_Ptl_BC1_type==1) then ! delete out-of_boundary particles from the list

            !$acc parallel loop
            do m=m0,m1
               zeon_a_arr(m)=zeon(1,m)
            enddo

            lost_num=0
            e_lost_num_BDY0=0
            e_lost_num_BDY1=0
            ! Find particle indices exiting the boundary
            do m=m0,m1
               if(Flag_Ptl_BC1_type==1 .and. zeon_a_arr(m) > psimax)then
                  lost_num = lost_num +1
                  lost_mid_arr(lost_num) = m
                  e_lost_num_BDY1 = e_lost_num_BDY1 +1
               elseif(Flag_Ptl_BC0_type==1 .and. zeon_a_arr(m) < psimin)then
                  lost_num = lost_num +1
                  lost_mid_arr(lost_num) = m
                  e_lost_num_BDY0 = e_lost_num_BDY0 +1
               endif
            enddo


            write(stdout,'(5(A,I8))') "istep=",istep,"; mype=",mype,"; e_lost_num=",lost_num, "=",e_lost_num_BDY0, "+",e_lost_num_BDY1
            ! write(stdout,*) "istep=",istep,"; mype=",mype,"electron lost_mid_arr=",lost_mid_arr(1:lost_num)
            call flush(stdout)


            if(lost_num>0) then
               call quicksort_integer(lost_mid_arr,1,lost_num)
            endif

            ! NOTE : mtop always should be mi even if m1!=mi.
            ! fill holes
            mtop = me
            me = me - lost_num ! remained me after loss
            lasth = lost_num

            allocate(fill_mid_arr(lost_num))
            fill_num=0
            do i=1, lost_num
               m = lost_mid_arr(i)
               if(m>me) exit
               do while(mtop == lost_mid_arr(lasth))
                  mtop=mtop-1
                  lasth=lasth-1
               enddo
               fill_mid_arr(i) = mtop
               mtop = mtop -1
               fill_num=i
            enddo

            ! write(stdout,*) "mype=",mype,"fill_num=",fill_num, "fill_mid_arr=",fill_mid_arr(1:fill_num)
            ! call flush(stdout)

            !$acc parallel loop
            do i=1, fill_num
               m = lost_mid_arr(i)
               mtop= fill_mid_arr(i)
               !$acc loop seq
               do j=1,dim_zeon
                  zeon(j, m)   = zeon(j, mtop)
               enddo
               !$acc loop seq
               do j=1,dim_zeon0
                  zeon0(j, m)  = zeon0(j, mtop)
               enddo
               elabels(m)      = elabels(mtop)
               we_init(m)      = we_init(mtop)
               we_vpara(m)     = we_vpara(mtop)
               we_vperp(m)     = we_vperp(mtop)
               we_accel(m)     = we_accel(mtop)
               we_adiab(m)     = we_adiab(mtop)
               we_fdamp(m)     = we_fdamp(mtop)

               UnNorm_theta_e(m) =UnNorm_theta_e(mtop)
               UnNorm_zeta_e(m) =UnNorm_zeta_e(mtop)
               if(Flag_RK4==1) then
                  !$acc loop seq
                  do j=1,dim_RK_cum_zeon
                     RK_cum_zeon(j,m)= RK_cum_zeon(j,mtop)
                  enddo
               endif
            enddo
            deallocate(fill_mid_arr)
            acc_newer_D_zeon=.true.

            call electron_locator_partial(m_start=m0, m_end=me)
            call n_T_u_group_e(m_start=m0, m_end=me)
         endif

      end subroutine treat_electrons_out_of_boundary


      subroutine treat_ions_out_of_boundary(m_start, m_end)

         use global_parameters
         use particle_array
         use field_array
         use function_gg_module
         use Ptl_locators_partial
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0, m1

         integer m,i,j
         real(wp) psimax,psimin
         integer lost_mid_arr(mimax) ! index of loss particles
         integer lost_mid_arr_BDY1(mimax)
         integer lost_num ! number of loss particles across boundary
         integer, allocatable :: fill_mid_arr(:)
         integer fill_num
         integer mtop, lasth
         integer dim_zion, dim_zion0, dim_RK_cum_zion
         real(wp) zion_a_arr(mi)
         real(wp), allocatable :: out_zion_BDY1(:,:), out_zion0_BDY1(:,:), out_ilabels_BDY1(:)

         m0=1
         m1=mi
         if(present(m_start)) m0=m_start
         if(present(m_end))   m1=m_end

         if(acc_newer_H_zion) then
            !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_H_zion=.false.
            acc_newer_D_zion=.false.
         endif

         dim_zion = size(zion,1)
         dim_zion0 = size(zion0,1)
         dim_RK_cum_zion = size(RK_cum_zion,1)

         psimax=psi_grid(mpsi)
         psimin=psi_grid(0)

         if(Flag_Ptl_BC0_type==0 .or. Flag_Ptl_BC1_type==0) then
            !$acc parallel loop present(zion,zion0)
            do m=m0,m1
               if(Flag_Ptl_BC1_type==0 .and. zion(1,m) > psimax)then
                  zion(1,m)=psimax
                  zion(2,m)=2.0_wp*pi-zion0(2,m)
                  zion(3,m)=zion0(3,m)
                  zion(4,m)=zion0(4,m)
                  zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zion(5,m)=-zion0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
                  zion(7,m)=zion0(7,m)
               elseif(Flag_Ptl_BC0_type==0 .and. zion(1,m) < psimin)then
                  zion(1,m)=psimin
                  zion(2,m)=2.0_wp*pi-zion0(2,m)
                  zion(3,m)=zion0(3,m)
                  zion(4,m)=zion0(4,m)
                  zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zion(5,m)=-zion0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
                  zion(7,m)=zion0(7,m)
               endif
            enddo
            acc_newer_D_zion=.true.

         endif

         if(Flag_Ptl_BC0_type==1 .or. Flag_Ptl_BC1_type==1) then ! delete out-of_boundary particles from the list

            !$acc parallel loop
            do m=m0,m1
               zion_a_arr(m)=zion(1,m)
            enddo

            lost_num=0
            i_lost_num_BDY0=0
            i_lost_num_BDY1=0
            ! Find particle indices exiting the boundary
            do m=m0,m1
               if(Flag_Ptl_BC1_type==1 .and. zion_a_arr(m) > psimax)then
                  lost_num = lost_num +1
                  lost_mid_arr(lost_num) = m
                  i_lost_num_BDY1 = i_lost_num_BDY1 +1
                  lost_mid_arr_BDY1(i_lost_num_BDY1) = m
               elseif(Flag_Ptl_BC0_type==1 .and. zion_a_arr(m) < psimin)then
                  lost_num = lost_num +1
                  lost_mid_arr(lost_num) = m
                  i_lost_num_BDY0 = i_lost_num_BDY0 +1
               endif
            enddo

            if(Flag_Recycling==1) then
               if(i_lost_num_BDY1>0) then
                  allocate(out_zion_BDY1(ziondim,i_lost_num_BDY1),out_zion0_BDY1(ziondim,i_lost_num_BDY1),out_ilabels_BDY1(i_lost_num_BDY1))

                  !$acc parallel loop copyout(out_zion_BDY1(:,:i_lost_num_BDY1),out_zion0_BDY1(:,:i_lost_num_BDY1),out_ilabels_BDY1(:i_lost_num_BDY1))
                  do m=1,i_lost_num_BDY1
                     out_zion_BDY1(:,m) = zion(:,lost_mid_arr_BDY1(m))
                     out_zion0_BDY1(:,m) = zion0(:,lost_mid_arr_BDY1(m))
                     out_ilabels_BDY1(m) = ilabels(lost_mid_arr_BDY1(m))
                  enddo

                  do m=1,i_lost_num_BDY1
                     if(out_zion0_BDY1(8,m)>0) then ! Ignore zero ftot particles
                        outION%mout = outION%mout +1

                        outION%P(:,outION%mout) = out_zion_BDY1(:,m)
                        outION%P0(:,outION%mout) = out_zion0_BDY1(:,m)
                        outION%labels(outION%mout) = out_ilabels_BDY1(m)
                     endif
                  enddo

                  deallocate(out_zion_BDY1,out_zion0_BDY1,out_ilabels_BDY1)
               endif
            endif


            write(stdout,'(5(A,I8))') "istep=",istep,"; mype=",mype,"; i_lost_num=",lost_num, "=",i_lost_num_BDY0, "+",i_lost_num_BDY1
            ! write(stdout,*) "istep=",istep,"; mype=",mype,"electron lost_mid_arr=",lost_mid_arr(1:lost_num)
            call flush(stdout)

            if(lost_num>0) then
               call quicksort_integer(lost_mid_arr,1,lost_num)
            endif

            ! NOTE : mtop always should be mi even if m1!=mi.
            ! fill holes
            mtop = mi
            mi = mi - lost_num ! remained mi after loss
            lasth = lost_num

            allocate(fill_mid_arr(lost_num))
            fill_num=0
            do i=1, lost_num
               m = lost_mid_arr(i)
               if(m>mi) exit
               do while(mtop == lost_mid_arr(lasth))
                  mtop=mtop-1
                  lasth=lasth-1
               enddo
               fill_mid_arr(i) = mtop
               mtop = mtop -1
               fill_num=i
            enddo

            ! write(stdout,*) "mype=",mype,"fill_num=",fill_num, "fill_mid_arr=",fill_mid_arr(1:fill_num)
            ! call flush(stdout)

            !$acc parallel loop
            do i=1, fill_num
               m = lost_mid_arr(i)
               mtop= fill_mid_arr(i)
               !$acc loop seq
               do j=1,dim_zion
                  zion(j, m)   = zion(j, mtop)
               enddo
               !$acc loop seq
               do j=1,dim_zion0
                  zion0(j, m)  = zion0(j, mtop)
               enddo
               ilabels(m)      = ilabels(mtop)
               wi_init(m)      = wi_init(mtop)
               wi_vpara(m)     = wi_vpara(mtop)
               wi_vperp(m)     = wi_vperp(mtop)
               wi_accel(m)     = wi_accel(mtop)
               wi_fdamp(m)     = wi_fdamp(mtop)

               UnNorm_theta_i(m) =UnNorm_theta_i(mtop)
               UnNorm_zeta_i(m) =UnNorm_zeta_i(mtop)
               if(Flag_RK4==1) then
                  !$acc loop seq
                  do j=1,dim_RK_cum_zion
                     RK_cum_zion(j,m)= RK_cum_zion(j,mtop)
                  enddo
               endif
            enddo
            deallocate(fill_mid_arr)
            acc_newer_D_zion=.true.

            call ion_locator_at_larmor_points_partial(m_start=m0, m_end=mi)
            call n_T_u_group(m_start=m0, m_end=mi)
         endif

      end subroutine treat_ions_out_of_boundary


      subroutine Pick_Random_Between_Numbers_NoDuplicates(a,b,Npick,out_arr)
         implicit none
#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"
         integer, intent(in) :: a,b,Npick ! Pick N random numbers between a and b without duplicates
         integer, intent(out) :: out_arr(:)
         integer i, r, len
         integer, allocatable :: nums_in_bin(:)

         len=b-a+1
         allocate(nums_in_bin(len))
         do i=1,len
            nums_in_bin(i)=a+(i-1) ! target numbers [a,b] in ascending order
         enddo

         do i=1,Npick
            r = floor(sprng()*len)+1 ! pick random number between [1,len]
            out_arr(i)=nums_in_bin(r)
            nums_in_bin(r) = nums_in_bin(len)
            len=len-1
         enddo

         deallocate(nums_in_bin)
      end subroutine Pick_Random_Between_Numbers_NoDuplicates



      function Produce_Recycling_Ion_and_Electrons_using_OutParticles() result(num_of_recycling_ptls)
         use global_parameters
         use particle_array
         use function_gg_module
         implicit none
#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"
         integer, intent(out) :: num_of_recycling_ptls
         integer m,k,ke,ki
         integer me0,mi0
         real(wp) :: ErgDampRatio=1e-2
         integer, allocatable :: rand_arr(:)


         mi0=mi+1 ! starting index of newly produced particles

         do k=1,outION%mout
            if(sprng()<outION%Recycling_Ratio) then
               mi=mi+1

               ! Copy to zion
               zion(:,mi)  = outION%P(:,k)
               zion0(:,mi) = outION%P0(:,k)
               ilabels(mi) = 1.0

               ! Reverse the parallel velocity and damp the magnitude
               zion(4,mi) = -sqrt(ErgDampRatio)*zion(4,mi)
               zion(6,mi) = sqrt(ErgDampRatio)*zion(6,mi)
               ! Reset the particle position to be the previous one
               zion(1:3,mi) = zion0(1:3,mi)
               zion(7,mi) = zion0(7,mi)

               ! For test
               ! zion(1,mi) = 0.5
               ! zion0(1,mi) = 0.5
               ! zion(7,mi) = 0.5
               ! zion0(7,mi) = 0.5
            endif
         enddo
         num_of_recycling_ptls = mi-mi0+1
         !$acc update device(zion(:,mi0:mi),zion0(:,mi0:mi),ilabels(mi0:mi))

         call n_T_u_group(mi0,mi,flag_backf0=.true.)
         call ion_locator_at_larmor_points_partial(mi0,mi)

         ! Determine the weight of the new particle
         !$acc parallel loop present(zion,zion0,wi_init,wi_vpara,wi_vperp,wi_accel,wi_fdamp,UnNorm_theta_i,UnNorm_zeta_i)
         do m=mi0,mi
            zion(5,m) = zion0(8,m)
            zion0(6,m) = 1e-6*zion0(6,m) ! Reduce the phase volume to reduce the particle noise
            zion0(8,m) = zion(5,m) + back_f0_i(m)*zion0(6,m)
         enddo

         !$acc parallel loop present(zion,zion0,wi_init,wi_vpara,wi_vperp,wi_accel,wi_fdamp,UnNorm_theta_i,UnNorm_zeta_i)
         do m=mi0,mi
            wi_init(m)  = zion(5,m)
            wi_vpara(m) = 0
            wi_vperp(m) = 0
            wi_accel(m) = 0
            wi_fdamp(m) = 0

            zion0(1:5,m)     = zion(1:5,m)
            zion0(7,m)       = zion(7,m)
            UnNorm_theta_i(m)= zion(2,m)
            UnNorm_zeta_i(m) = zion(3,m)
         enddo



         ! Produce electrons using the new ion particles
         me0=me+1
         me =me+num_of_recycling_ptls


         !$acc parallel loop present(zeon,zion,backTe,backT)
         do k=1,num_of_recycling_ptls
            ke=me0+(k-1)
            ki=mi0+(k-1)

            zeon(:,ke) = zion(:,ki)
            zeon0(:,ke) = zion0(:,ki)

            ! Note that zion0(8,m) is the full-f value, but zion(8,m) is just a temporary old_backf0 value
            ! whereas zeon(8,m) and zeon0(8,m) are the full-h (or full-f) marker values at different time step
            ! Therefore, we must set zeon(8,m) to be the same as the zeon0(8,m) (not as the same as zion(8,m) which is trash value here)
            ! since zeon0(8,m) will be overwritten by zeon(8,m) at the push routines
            zeon(8,ke) = zeon0(8,ke)
            ! elabels(ke) = ilabels(ki) ! NOTE: This is not correct.. just for test now
            elabels(ke) = 1.0
         enddo

         call n_T_u_group_e(me0,me,flag_backf0=.true.)  ! get Te, note that back_f0_e cannot be calculated here because it requires zion(4) and zion(6)


         allocate(rand_arr(num_of_recycling_ptls))
         do k=1,num_of_recycling_ptls
            rand_arr(k) = sprng()
         enddo

         !$acc parallel loop present(zeon,zion,backTe,backT)
         do k=1,num_of_recycling_ptls
            ke=me0+(k-1)
            ki=mi0+(k-1)

            zeon(4,ke)=zion(4,ki)*sqrt(msm0)*c0cs*sqrt(msm0*backTe(ke)/backT(ki))  !checked; neglect electron rotation
            zeon(6,ke)=zion(6,ki)*sqrt(msm0*backTe(ke)/backT(ki))     ! =sqrt(mu_e)

            ! Electron's vpara direction is random because of the ionization process
            zeon(4,ke)=sign(1.0,rand_arr(k)-0.5)*zeon(4,ke)
         enddo

         call electron_locator_partial(me0,me) ! This is needed to calculate evector (gradPhi) and its interpolation (wpi)

         !$acc parallel loop present(zeon,zeon0,we_init,we_vpara,we_vperp,we_accel,we_adiab,we_fdamp,UnNorm_theta_e,UnNorm_zeta_e)
         do m=me0,me
            we_init(m)  = zeon(5,m)
            we_vpara(m) = 0
            we_vperp(m) = 0
            we_accel(m) = 0
            we_adiab(m) = 0
            we_fdamp(m) = 0

            zeon0(1:5,m)     = zeon(1:5,m)
            zeon0(7,m)       = zeon(7,m)
            UnNorm_theta_e(m)= zeon(2,m)
            UnNorm_zeta_e(m) = zeon(3,m)
         enddo


         call Clear_OutParticles

         deallocate(rand_arr)

      end function Produce_Recycling_Ion_and_Electrons_using_OutParticles


end module Ptl_Src_Loss
