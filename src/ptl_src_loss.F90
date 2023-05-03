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


   type Particle_Source_Control_Variables
      integer :: dim_a, dim_gq
      real(wp), allocatable, dimension(:) :: a_1D, gq_1D, s_1D, Te0_1D, Te0_eV_1D, Ti0_1D, Ti0_eV_1D, safety_q_1D
      real(wp), allocatable, dimension(:,:) :: Cell_dV
      real(wp) :: a_min, a_max, da, dgq, inv_da, inv_dgq

      real(wp) :: ptl_label ! laberl of particles generated within each PSCV
      integer :: nstep, cum_num
      real(wp), allocatable, dimension(:,:,:) :: Elec_dupara_per_nstep,Ion_dupara_per_nstep
      real(wp), allocatable, dimension(:,:,:) :: Elec_dupara_per_nstep_SM,Ion_dupara_per_nstep_SM

      real(wp), allocatable, dimension(:,:,:) :: Elec_dErg_eV_per_nstep, Ion_dErg_eV_per_nstep
   end type

   type(OutParticles) outION
   type(Particle_Source_Control_Variables), pointer :: PSCV
   type(Particle_Source_Control_Variables), target :: PSCV_list(3) ! (1:Inside,  2: beyond BDY0, 3: beyond BDY1)

   contains

      subroutine Initialize_Ptl_Src_Loss
         use global_parameters
         use Structured_Grid
         implicit none

         if(Flag_Recycling==1) then
            outION%mout_max = mimax/20
            allocate(outION%P(ziondim,outION%mout_max),outION%P0(ziondim,outION%mout_max))
            allocate(outION%energy(outION%mout_max),outION%backT(outION%mout_max),outION%backf0(outION%mout_max),outION%labels(outION%mout_max))
         endif

         call Initialize_FootPrint(FP_BDY1,a1)
         call Initialize_FootPrint(cum_FP_BDY1,a1)

         PSCV=>PSCV_list(1) ! beyond inner boundary (BDY0)
         call Initialize_Particle_Source_Control_Variables(minval(SG%a_1D)-0.05, minval(SG%a_1D), 10, Ptl_Src_dim_gq, -1.0)
         PSCV=>PSCV_list(2) ! Inside the domain
         call Initialize_Particle_Source_Control_Variables(minval(SG%a_1D), maxval(SG%a_1D), Ptl_Src_dim_a, Ptl_Src_dim_gq, -2.0)
         PSCV=>PSCV_list(3) ! beyond outer boundary (BDY1)
         call Initialize_Particle_Source_Control_Variables(maxval(SG%a_1D), maxval(SG%a_1D)+0.07, 20, Ptl_Src_dim_gq, -3.0)

      end subroutine Initialize_Ptl_Src_Loss

      subroutine Initialize_Particle_Source_Control_Variables(in_a_min, in_a_max, in_dim_a, in_dim_gq, in_ptl_label)
         use global_parameters
         use Structured_Grid
         use function_gg_module
         implicit none
         real(wp),intent(in) :: in_a_min, in_a_max
         integer, intent(in) :: in_dim_a, in_dim_gq
         real(wp),intent(in) :: in_ptl_label
         integer :: dim_gq, dim_a, aid, gqid
         real(wp) :: this_a,this_da, this_gq
         real(wp), external :: Jacobian

         PSCV%a_min = in_a_min
         PSCV%a_max = in_a_max
         PSCV%dim_a = in_dim_a
         PSCV%dim_gq= in_dim_gq
         PSCV%ptl_label = in_ptl_label

         dim_a=PSCV%dim_a
         dim_gq=PSCV%dim_gq
         allocate(PSCV%a_1D(dim_a),PSCV%gq_1D(dim_gq), PSCV%s_1D(dim_a))
         allocate(PSCV%Te0_1D(dim_a), PSCV%Te0_eV_1D(dim_a), PSCV%Ti0_1D(dim_a), PSCV%Ti0_eV_1D(dim_a), PSCV%safety_q_1D(dim_a))
         allocate(PSCV%Cell_dV(dim_gq,dim_a))

         do aid=1,PSCV%dim_a
            PSCV%a_1D(aid) = PSCV%a_min  +(aid-1)*(PSCV%a_max-PSCV%a_min)/(PSCV%dim_a-1)
         enddo
         do gqid=1,PSCV%dim_gq
            PSCV%gq_1D(gqid) = (gqid-1)*(2.0*pi)/(PSCV%dim_gq-1)
         enddo
         call Convert_a_arr_to_s_arr(PSCV%a_1D, PSCV%s_1D)

         call Interpolate_SG_1D_Data_at_given_s_arr(SG%Te0_1D, PSCV%s_1D, PSCV%Te0_1D)
         call Interpolate_SG_1D_Data_at_given_s_arr(SG%Ti0_1D, PSCV%s_1D, PSCV%Ti0_1D)
         call Interpolate_SG_1D_Data_at_given_s_arr(SG%Te0_eV_1D, PSCV%s_1D, PSCV%Te0_eV_1D)
         call Interpolate_SG_1D_Data_at_given_s_arr(SG%Ti0_eV_1D, PSCV%s_1D, PSCV%Ti0_eV_1D)
         call Interpolate_SG_1D_Data_at_given_s_arr(SG%safety_q_1D, PSCV%s_1D, PSCV%safety_q_1D)


         PSCV%da = PSCV%a_1D(2)-PSCV%a_1D(1)
         PSCV%dgq = PSCV%gq_1D(2)-PSCV%gq_1D(1)

         PSCV%inv_da = 1.0/PSCV%da
         PSCV%inv_dgq = 1.0/PSCV%dgq
         do aid=1,PSCV%dim_a
            do gqid=1,PSCV%dim_gq
               this_gq= PSCV%gq_1D(gqid)

               this_da = min(PSCV%a_max,(PSCV%a_1D(aid)+0.5*PSCV%da)) - max(PSCV%a_min,PSCV%a_1D(aid)-0.5*PSCV%da)
               this_a  = 0.5*(min(PSCV%a_max,(PSCV%a_1D(aid)+0.5*PSCV%da))+max(PSCV%a_min,PSCV%a_1D(aid)-0.5*PSCV%da))

               PSCV%Cell_dV(gqid,aid) = Jacobian(this_a,this_gq)*this_da*PSCV%dgq*(zetamax-zetamin)
            enddo
         enddo


         PSCV%nstep = Ptl_Src_nstep
         allocate(PSCV%Elec_dupara_per_nstep(dim_gq,dim_a,0:1), PSCV%Ion_dupara_per_nstep(dim_gq,dim_a,0:1))
         allocate(PSCV%Elec_dupara_per_nstep_SM(dim_gq,dim_a,0:1), PSCV%Ion_dupara_per_nstep_SM(dim_gq,dim_a,0:1))
         allocate(PSCV%Elec_dErg_eV_per_nstep(dim_gq,dim_a,0:1),PSCV%Ion_dErg_eV_per_nstep(dim_gq,dim_a,0:1))

         PSCV%cum_num = 0
      end subroutine Initialize_Particle_Source_Control_Variables

      subroutine Finalize_Ptl_Src_Loss
         implicit none
         integer i

         if(allocated(outION%P)) deallocate(outION%P, outION%P0, outION%energy, outION%backT, outION%backf0,outION%labels)

         do i=1,size(PSCV_list)
            PSCV=>PSCV_list(i)
            if(allocated(PSCV%Elec_dupara_per_nstep)) deallocate(PSCV%Elec_dupara_per_nstep,PSCV%Ion_dupara_per_nstep)
            if(allocated(PSCV%Elec_dupara_per_nstep_SM)) deallocate(PSCV%Elec_dupara_per_nstep_SM,PSCV%Ion_dupara_per_nstep_SM)
            if(allocated(PSCV%Elec_dErg_eV_per_nstep)) deallocate(PSCV%Elec_dErg_eV_per_nstep,PSCV%Ion_dErg_eV_per_nstep)
            if(allocated(PSCV%a_1D)) deallocate(PSCV%a_1D,PSCV%gq_1D,PSCV%s_1D)
            if(allocated(PSCV%Te0_eV_1D)) deallocate(PSCV%Te0_1D, PSCV%Te0_eV_1D, PSCV%Ti0_1D, PSCV%Ti0_eV_1D, PSCV%safety_q_1D)
            if(allocated(PSCV%Cell_dV)) deallocate(PSCV%Cell_dV)
         enddo

      end subroutine Finalize_Ptl_Src_Loss

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


      subroutine Cal_dupara_per_dt_using_SG
         use global_parameters
         use Structured_Grid
         implicit none
         real(wp) dt_s,me_kg,mi_kg,qe_C,qi_C, vthe_mps, vthi_mps ! velocity unit (meter per second)
         integer aa,qq,k
         real(wp) peak_before, peak_after
         integer Wsz_q, Wsz_a

         dt_s  = tstep*utime
         me_kg = 1e-3*amasse
         mi_kg = 1e-3*amassi
         qe_C = charge/2997924580.0
         qi_C = chargi/2997924580.0

         if(istep>PSCV%nstep .and. mod(istep,PSCV%nstep)/=0) then
            write(stdout,'(2(A,I0))') "Something is Wrong!! istep=", istep,",  is not multiple of PtlSrcCV%nstep=", PSCV%nstep
            call flush(stdout)
            stop
         else

            do aa=1,PSCV%dim_a
               vthe_mps = sqrt(abs(qe_C)*PSCV%Te0_eV_1D(aa)/(me_kg))
               vthi_mps = sqrt(abs(qi_C)*PSCV%Ti0_eV_1D(aa)/(mi_kg))
            enddo
         endif

         Wsz_q=3;
         Wsz_a=2;
         do aa=1,PSCV%dim_a
            do qq=1,PSCV%dim_gq
               do k=0,1
                  PSCV%Elec_dupara_per_nstep_SM(qq,aa,k) = maxval(PSCV%Elec_dupara_per_nstep(modulo([qq-Wsz_q:qq+Wsz_q]-1, PSCV%dim_gq)+1, max(1,aa-Wsz_a):min(PSCV%dim_a,aa+Wsz_a),k))
                  PSCV%Ion_dupara_per_nstep_SM(qq,aa,k) = maxval(PSCV%Ion_dupara_per_nstep(modulo([qq-Wsz_q:qq+Wsz_q]-1, PSCV%dim_gq)+1, max(1,aa-Wsz_a):min(PSCV%dim_a,aa+Wsz_a),k))
               enddo
            enddo
         enddo

      end subroutine Cal_dupara_per_dt_using_SG


      subroutine Cal_dErg_per_nstep_by_Potential_along_FieldLine
         use global_parameters
         use Magnetic_Grid, only: B_cntv,grad_B_cntv,d23_B_cntv,TriCubic_Interpolation_of_Vec, MG_a_1D=>a_1D
         use particle_decomp
         use mpi_defs
         use Structured_Grid
         use function_gg_module
         implicit none

         integer :: minpts, ifail=1, ierror, istatus(MPI_STATUS_SIZE)
         real(wp) :: acc, error_8=1e-5
         real(wp), external:: dsda, D01AHF
         real(wp) dgf,specific_time ! toroidal step for field-line-following integration
         real(wp) a_RK(4), gq_RK(4), B_vec_cntv(3)
         integer, allocatable :: F_ini_aid(:), F_ini_gqid(:)
         real(wp),dimension(:),allocatable :: ini_a,ini_gq,ini_gf
         real(wp) now_a, now_gq, now_gf, now_s, ori_phi_V, now_phi_V ! Note that phi is the potential in Volt unit
         real(wp) prev_a, prev_gq, prev_gf
         real(wp):: aaa, the, Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj
         integer recv_counts_Lsize(npartdom), recv_displs_Lsize(npartdom)
         integer max_tor_rotation, step,max_step,num_gq, num_gf, i, j, k, ierr,length, h5err
         integer Gsize, Lsize, extra! Number of G(Global) or L(Local) Field lines
         integer Fid, Fid_start, Fid_end
         integer tor_dir,aid,gqid
         real(wp) :: dErg_eV, FL_loc_Lc ! field line length
         real(wp), allocatable :: FL_max_Elec_dErg_eV(:),FL_max_Ion_dErg_eV(:), recv_FL_max_dErg_eV(:) ! maximum energy change by potential along field line
         real(wp) Lc_track_threshold ! Determine how long the field line will be tracked (vmax*dt*nstep)
         real(wp) Elec_Lc_thold, Ion_Lc_thold
         real(wp) MG_a_min, MG_a_max

#ifdef _MGY_TIMER_
   Gt_start=MPI_WTIME()
#endif

         MG_a_min = minval(MG_a_1D)
         MG_a_max = maxval(MG_a_1D)

         ! Prepare SG's phi_3D in case
         call Update_SG_phi_3D_by_gathering_SG_phi2D

         specific_time = 0.0


         max_tor_rotation = 100
         ! Lc_track_threshold = 50 ! 50 m

         dgf = 2.0*pi/100.0
         Elec_Lc_thold = PSCV%nstep*(2.0*umax*MM_vthce*tstep) ! electron travel distance (m) during nstep assuming 2*umax speed
         Ion_Lc_thold  = PSCV%nstep*(2.0*umax*MM_vthci*tstep) ! ion travel distance (m) during nstep assuming 2*umax speed
         Lc_track_threshold = max(Elec_Lc_thold, Ion_Lc_thold)
         max_step = ceiling(Lc_track_threshold/(MM_R_axis*dgf)) ! Changle Lc condition into approximate number of steps


         if(mype==0) then
            write(stdout,'(A,I0,2(A,ES8.2),A,I0)') "mype=",mype,"; Elec_Lc_thold =", Elec_Lc_thold, " (m), Ion_Lc_thold=", Ion_Lc_thold, " (m); max_step=",max_step
            call flush(stdout)
         endif


         Gsize = PSCV%dim_a*PSCV%dim_gq
         Lsize = Gsize/npartdom
         extra = Gsize-Lsize*npartdom
         Fid_start = Lsize*(myrank_partd)   + min(extra,myrank_partd) +1 ! Last +1 is due to Fortran's indexing that starts from 1
         Fid_end   = Lsize*(myrank_partd+1) + min(extra,myrank_partd+1)
         Lsize = Fid_end - Fid_start +1

         ! Get counts and displacements of Lsize in MPI_COMM_WORLD
         call MPI_ALLGATHER(Lsize,1,MPI_INTEGER,recv_counts_Lsize,1,MPI_INTEGER,partd_comm,ierr)
         recv_displs_Lsize(1)=0
         do i=1,npartdom-1
            recv_displs_Lsize(i+1)=recv_displs_Lsize(i)+recv_counts_Lsize(i)
         enddo


         allocate(F_ini_aid(Gsize), F_ini_gqid(Gsize))
         allocate(ini_a(Gsize),ini_gq(Gsize),ini_gf(Gsize))

         ! Calculate at zetamax and communicate with other planes later
         ini_gf(:) = zetamax
         Fid = 1
         do aid=1,PSCV%dim_a
            do gqid=1,PSCV%dim_gq
               ini_a(Fid) = PSCV%a_1D(aid)
               ini_gq(Fid)= PSCV%gq_1D(gqid)
               F_ini_aid(Fid) = aid
               F_ini_gqid(Fid)= gqid
               Fid=Fid+1
            enddo
         enddo


         allocate(FL_max_Elec_dErg_eV(Fid_start:Fid_end), FL_max_Ion_dErg_eV(Fid_start:Fid_end))
         allocate(recv_FL_max_dErg_eV(Gsize))

         FL_max_Elec_dErg_eV = 0.0
         FL_max_Ion_dErg_eV  = 0.0
         !$acc parallel loop private(a_RK,gq_RK,B_vec_cntv)
         do Fid = Fid_start,Fid_end
            !$acc loop seq private(a_RK,gq_RK)
            do tor_dir=-1,1,2

               FL_loc_Lc = 0.0

               ! dgf = tor_dir*2.0*pi/50.0 ! tor_dir=-1 (CW; -zeta), tor_dir=1(CCW; +zeta)
               dgf = tor_dir*abs(dgf) ! tor_dir=-1 (CW; -zeta), tor_dir=1(CCW; +zeta)

               now_a  = ini_a(Fid)
               now_gq = ini_gq(Fid)
               now_gf = ini_gf(Fid)

               prev_a  = now_a
               prev_gq = now_gq
               prev_gf = now_gf
               call Convert_a_to_s(now_a,now_s)
               call Get_Intp_global_phi3D_on_SG(now_s, now_gq, now_gf, ori_phi_V)

               ! do while(.true.)
               !$acc loop seq private(a_RK,gq_RK)
               do step=1,max_step

                  prev_a = now_a
                  prev_gq = now_gq
                  prev_gf = now_gf

                  ! Runge-Kutta 4th order
                  ! call esiget_Intp_B_vec_contravariant(now_a, modulo(now_gq,2.0*pi), modulo(now_gf,2.0*pi), specific_time, B_vec_cntv)
                  call TriCubic_Interpolation_of_Vec(now_a, modulo(now_gq,2.0*pi), modulo(now_gf,2.0*pi), B_cntv, grad_B_cntv, d23_B_cntv, B_vec_cntv)
                  a_RK(1)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
                  gq_RK(1)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
                  ! call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(1), modulo(now_gq+0.5*gq_RK(1),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), specific_time, B_vec_cntv)
                  call TriCubic_Interpolation_of_Vec(now_a+0.5*a_RK(1), modulo(now_gq+0.5*gq_RK(1),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), B_cntv, grad_B_cntv, d23_B_cntv, B_vec_cntv)
                  a_RK(2)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
                  gq_RK(2)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
                  ! call esiget_Intp_B_vec_contravariant(now_a+0.5*a_RK(2), modulo(now_gq+0.5*gq_RK(2),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), specific_time, B_vec_cntv)
                  call TriCubic_Interpolation_of_Vec(now_a+0.5*a_RK(2), modulo(now_gq+0.5*gq_RK(2),2.0*pi), modulo(now_gf+0.5*dgf,2.0*pi), B_cntv, grad_B_cntv, d23_B_cntv, B_vec_cntv)
                  a_RK(3)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
                  gq_RK(3)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv
                  ! call esiget_Intp_B_vec_contravariant(now_a+a_RK(3), modulo(now_gq+gq_RK(3),2.0*pi), modulo(now_gf+dgf,2.0*pi), specific_time, B_vec_cntv)
                  call TriCubic_Interpolation_of_Vec(now_a+a_RK(3), modulo(now_gq+gq_RK(3),2.0*pi), modulo(now_gf+dgf,2.0*pi), B_cntv, grad_B_cntv, d23_B_cntv, B_vec_cntv)
                  a_RK(4)     = dgf*B_vec_cntv(1)/B_vec_cntv(3) ! Ba_cntv/Bgf_cntv
                  gq_RK(4)    = dgf*B_vec_cntv(2)/B_vec_cntv(3) ! Bgq_cntv/Bgf_cntv


                  now_a = now_a + (a_RK(1)+2.0*a_RK(2)+2.0*a_RK(3)+a_RK(4))/6.0
                  now_gq = now_gq + (gq_RK(1)+2.0*gq_RK(2)+2.0*gq_RK(3)+gq_RK(4))/6.0
                  now_gf = now_gf + dgf


                  ! aaa=now_a
                  ! the=modulo(now_gq, 2.0*pi)
                  ! call esiget2dthese(Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj,aaa,the,1)
                  ! call esiget_Intp_B_vec_contravariant(aaa, the, modulo(now_gf,2.0*pi), specific_time, B_vec_cntv)
                  ! call TriCubic_Interpolation_of_Vec(aaa, the, modulo(now_gf,2.0*pi), B_cntv, grad_B_cntv, d23_B_cntv, B_vec_cntv)
                  ! FL_loc_Lc  = FL_loc_Lc + abs(dgf)*sqrt( ((B_vec_cntv(1)/B_vec_cntv(3))**2)*(raj**2+zaj**2)  &
                  !                                                    +((B_vec_cntv(2)/B_vec_cntv(3))**2)*(rqj**2+zqj**2)  &
                  !                                                    +(B_vec_cntv(1)*B_vec_cntv(2)/(B_vec_cntv(3)**2))*(raj*rqj+zaj*zqj) &
                  !                                                    + rj**2  )


                  FL_loc_Lc = FL_loc_Lc+MM_R_axis*abs(dgf)

                  ! Interpolate phi at current point (now_a,now_gq,now_gf)

                  call Convert_a_to_s(now_a,now_s)
                  call Get_Intp_global_phi3D_on_SG(now_s, now_gq, now_gf, now_phi_V)

                  dErg_eV = abs(now_phi_V-ori_phi_V)
                  if(dErg_eV > FL_max_Elec_dErg_eV(Fid)) then
                     FL_max_Elec_dErg_eV(Fid) = dErg_eV
                  endif

                  if(FL_loc_Lc < Ion_Lc_thold .and. dErg_eV > FL_max_Ion_dErg_eV(Fid)) then
                     FL_max_Ion_dErg_eV(Fid) = dErg_eV
                  endif

                  ! if(now_a<=PSCV%a_min .or. now_a>=PSCV%a_max) then
                  if(now_a<=MG_a_min+1e-2 .or. now_a>=MG_a_max-1e-2) then
                     exit
                  endif
               enddo ! end while(.true.)
            enddo ! tor_id=1,2

            ! write(stdout,'(2(A,I0),A,2(I0,1x),2(A,ES9.3))') "mype=",mype,"; Fid=",Fid,"; Fid_start_end=",Fid_start, Fid_end, "; FL_max_dErg_eV(Fid)", FL_max_dErg_eV(Fid), "; FL_loc_Lc=", FL_loc_Lc
            ! call flush(stdout)
            ! call MPI_Barrier(MPI_COMM_WORLD,ierrr)
         enddo ! Fid

         ! Gather Electron max_dErg
         call MPI_Allgatherv(FL_max_Elec_dErg_eV, Lsize, MPI_Rsize, recv_FL_max_dErg_eV, recv_counts_Lsize, recv_displs_Lsize, MPI_Rsize, partd_comm, ierr)
         do Fid=1,Gsize
            PSCV%Elec_dErg_eV_per_nstep(F_ini_gqid(Fid),F_ini_aid(Fid),1) = recv_FL_max_dErg_eV(Fid)
         enddo
#ifdef _MGY_TIMER_
         if(mype==0) write(stdout,'(2(A,I0),2(A,ES11.3),A)') "I'm=",mype,", Ptl_Src_nstep=",Ptl_Src_nstep,", PSCV_mean_a=",0.5*(PSCV%a_min+PSCV%a_max),", max Elec_dErg_eV_per_nstep=", maxval(recv_FL_max_dErg_eV)," eV"
#endif

         ! Gather Ion max_dErg
         call MPI_Allgatherv(FL_max_Ion_dErg_eV, Lsize, MPI_Rsize, recv_FL_max_dErg_eV, recv_counts_Lsize, recv_displs_Lsize, MPI_Rsize, partd_comm, ierr)
         do Fid=1,Gsize
            PSCV%Ion_dErg_eV_per_nstep(F_ini_gqid(Fid),F_ini_aid(Fid),1) = recv_FL_max_dErg_eV(Fid)
         enddo
#ifdef _MGY_TIMER_
         if(mype==0) write(stdout,'(2(A,I0),2(A,ES11.3),A)') "I'm=",mype,", Ptl_Src_nstep=",Ptl_Src_nstep,", PSCV_mean_a=",0.5*(PSCV%a_min+PSCV%a_max),", max Ion_dErg_eV_per_nstep=", maxval(recv_FL_max_dErg_eV)," eV"
#endif
         !toroidal communication
         length = PSCV%dim_a*PSCV%dim_gq
         call MPI_SENDRECV(PSCV%Elec_dErg_eV_per_nstep(:,:,1),length,mpi_Rsize,left_pe,1, PSCV%Elec_dErg_eV_per_nstep(:,:,0),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         call MPI_SENDRECV(PSCV%Ion_dErg_eV_per_nstep(:,:,1),length,mpi_Rsize,left_pe,1, PSCV%Ion_dErg_eV_per_nstep(:,:,0),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)

         deallocate(F_ini_aid, F_ini_gqid)
         deallocate(ini_a, ini_gq, ini_gf)
         deallocate(FL_max_Elec_dErg_eV, FL_max_Ion_dErg_eV)
         deallocate(recv_FL_max_dErg_eV)

#ifdef _MGY_TIMER_
   Gt_end=MPI_WTIME()
   if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "I'm=",mype, ", Cal_dErg_per_nstep_by_Potential_along_FieldLine takes", Gt_end-Gt_start, " secs"
   call flush(stdout)
#endif

      end subroutine Cal_dErg_per_nstep_by_Potential_along_FieldLine


      subroutine Cal_dupara_per_dt_using_dErg_per_nstep
         use global_parameters
         implicit none
         real(wp) dt_s,me_kg,mi_kg,qe_C,qi_C, vthe_mps, vthi_mps, vmax_mps ! velocity unit (meter per second)
         real(wp) MarginFactor, MinimumVal
         integer aid,ss,qq,k
         integer Wsz_q, Wsz_a

         dt_s  = tstep*utime
         me_kg = 1e-3*amasse
         mi_kg = 1e-3*amassi
         qe_C = charge/2997924580.0
         qi_C = chargi/2997924580.0

         MarginFactor = 1.1 ! Give 10% margin for more safety
         MinimumVal = 0.05 ! du=0.05 everywhere

         do ss=1,PSCV%dim_a
            vthe_mps = sqrt(abs(qe_C)*PSCV%Te0_eV_1D(ss)/(me_kg))
            vthi_mps = sqrt(abs(qi_C)*PSCV%Ti0_eV_1D(ss)/(mi_kg))

            vmax_mps = umax*vthe_mps
            PSCV%Elec_dupara_per_nstep(:,ss,:) = max(MinimumVal, abs(-vmax_mps + sqrt(vmax_mps**2.0 + abs(qe_C)*MarginFactor*2.0*PSCV%Elec_dErg_eV_per_nstep(:,ss,:)/me_kg))/(vthe_mps))

            vmax_mps = umax*vthi_mps
            PSCV%Ion_dupara_per_nstep(:,ss,:)  = max(MinimumVal, abs(-vmax_mps + sqrt(vmax_mps**2.0 + abs(qe_C)*MarginFactor*2.0*PSCV%Ion_dErg_eV_per_nstep(:,ss,:)/mi_kg))/(vthi_mps))
         enddo

         ! Wsz_q=5;
         ! Wsz_a=3;
         Wsz_q=max(3,ceiling(0.02*PSCV%dim_gq))
         Wsz_a=max(2,ceiling(0.02*PSCV%dim_a))
         do ss=1,PSCV%dim_a
            do qq=1,PSCV%dim_gq
               do k=0,1
                  PSCV%Elec_dupara_per_nstep_SM(qq,ss,k) = max(MinimumVal, maxval(PSCV%Elec_dupara_per_nstep(modulo([qq-Wsz_q:qq+Wsz_q]-1, PSCV%dim_gq)+1, max(1,ss-Wsz_a):min(PSCV%dim_a,ss+Wsz_a),k)))
                  PSCV%Ion_dupara_per_nstep_SM(qq,ss,k) = max(MinimumVal, maxval(PSCV%Ion_dupara_per_nstep(modulo([qq-Wsz_q:qq+Wsz_q]-1, PSCV%dim_gq)+1, max(1,ss-Wsz_a):min(PSCV%dim_a,ss+Wsz_a),k)))
               enddo
            enddo
         enddo

      end subroutine Cal_dupara_per_dt_using_dErg_per_nstep


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
         if(Flag_Ptl_BC0_type==1) call load_electrons_In_or_Out_Flux(dNe_BDY0, delta_t_seconds,Ptl_BC0_val,a0)
         ! Create Src from inner boundary (BDY1)
         if(Flag_Ptl_BC1_type==1) call load_electrons_In_or_Out_Flux(dNe_BDY1, delta_t_seconds,Ptl_BC1_val,a1)
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
         if(Flag_Ptl_BC0_type==1) call load_ions_In_or_Out_Flux(dNi_BDY0, delta_t_seconds,Ptl_BC0_val,a0)
         ! Create Src from inner boundary (BDY1)
         if(Flag_Ptl_BC1_type==1) call load_ions_In_or_Out_Flux(dNi_BDY1, delta_t_seconds,Ptl_BC1_val,a1)
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

         integer :: Nm,m0,m,tmp_num,k
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

         real(wp), allocatable :: wpi(:,:),gaE(:),wgaE(:,:),gatE(:),dgr(:),da(:),dgq(:),dgf(:),Er_particle(:)
         !$acc declare create(dgr,da,dgq,dgf,wpi,Er_particle,gaE,wgaE,gatE)
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
         else
            write(stdout,*) "Error!! Unknown Flag; Flag_Load_and_WeightDef=",Flag_Load_and_WeightDef
            call flush(stdout)
            stop
         endif
         !$acc update device(zeon(:,m0:me))

         ! call n_T_u_group_e(m0)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))
         call n_T_u_group_e(m0,me,flag_backf0=.true.)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))
         !$acc update self(back_f0_e(m0:me))

         call electron_locator_partial(m0,me) ! This is needed to calculate evector (gradPhi) and its interpolation (wpi)





         allocate(wpi(3,me),gaE(me),wgaE(3,me),gatE(me),dgr(me),da(me),dgq(me),dgf(me),Er_particle(me))
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

         if(electromagnetic >= 1)then
            !$acc parallel loop present(zeon,gaE,wgaE,gatE)
            do m=m0,me
               call Get_MP_alpha_amplitude(zeon(1,m),zeon(2,m),zeon(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
            enddo

            if(Flag_GC_Solver==0) then
               m=gcmotion_eon_em_MGY(dgr,da,dgq,dgf,zeon,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,me,zeondim)
            elseif(Flag_GC_Solver==1) then
               m=gcmotion_eon_em_MGY_vec_form(dgr,da,dgq,dgf,zeon,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,m0,me,zeondim)
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
            ! zeon(8,m) is "full-f or full-h"
            ! zeon(9,m) is "f_AD"
            ! zeon0(9,m) is "f_tot"
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

            !$acc loop seq
            do k=1,5
               zeon0(k,m)    = zeon(k,m)
            enddo
            zeon0(7,m)       = zeon(7,m)
            zeon0(8,m)       = zeon(8,m)
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
         deallocate(wpi,gaE,wgaE,gatE,dgr,da,dgq,dgf,Er_particle)
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

         integer :: Nm,m0,m,tmp_num,k
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

         real(wp), allocatable :: wpi(:,:),gaE(:),wgaE(:,:),gatE(:),dgr(:),da(:),dgq(:),dgf(:),Er_particle(:)
         !$acc declare create(dgr,da,dgq,dgf,wpi,Er_particle,gaE,wgaE,gatE)
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
         else
            write(stdout,*) "Error!! Unknown Flag; Flag_Load_and_WeightDef=",Flag_Load_and_WeightDef
            call flush(stdout)
            stop
         endif
         !$acc update device(zion(:,m0:mi))

         ! call n_T_u_group(m0)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
         call n_T_u_group(m0,mi,flag_backf0=.true.)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
         !$acc update self(back_f0_i(m0:mi))


         call ion_locator_at_larmor_points_partial(m0,mi) ! This is needed to calculate evector (gradPhi) and its interpolation (wpi)



         allocate(wpi(3,mi),gaE(mi),wgaE(3,mi),gatE(mi),dgr(mi),da(mi),dgq(mi),dgf(mi),Er_particle(mi))



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

         if(electromagnetic >= 1)then
            !$acc parallel loop present(zion,gaE,wgaE,gatE)
            do m=m0,mi
               call Get_MP_alpha_amplitude(zion(1,m),zion(2,m),zion(3,m),istep*tstep, gaE(m), wgaE(:,m), gatE(m))
            enddo

            if(Flag_GC_Solver==0) then
               m=gcmotion_em_MGY(dgr,da,dgq,dgf,zion,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,mi,ziondim)
            elseif(Flag_GC_Solver==1) then
               m=gcmotion_em_MGY_vec_form(dgr,da,dgq,dgf,zion,wpi,gaE,wgaE,gatE,Er_particle,nonlinear,m0,mi,ziondim)
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

            !$acc loop seq
            do k=1,5
               zion0(k,m)     = zion(k,m)
            enddo
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
         deallocate(wpi,gaE,wgaE,gatE,dgr,da,dgq,dgf,Er_particle)

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
         integer m,i,j,k,is_this_trc
         real(wp) psimax,psimin
         integer lost_mid_arr(memax) ! index of loss particles
         integer lost_num ! number of loss particles across boundary
         integer, allocatable :: fill_mid_arr(:)
         integer fill_num
         integer mtop, lasth
         integer dim_zeon, dim_zeon0, dim_RK_cum_zeon
         real(wp) zeon_a_arr(me)
         real(wp) target_w_tot, integral_f0

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
                  zeon(1,m)=zeon0(1,m)
                  zeon(2,m)=2.0_wp*pi-zeon0(2,m)
                  ! Do not reset zeon(3,m) or toroidal angle !!
                  ! This is to ensure that this particle is in the same toroidal rank.
                  zeon(4,m)=zeon0(4,m)
                  zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
                  zeon(7,m)=zeon0(7,m)

                  if(Flag_w_Evol>=3) zeon(8,m)=zeon0(8,m)
               elseif(Flag_Ptl_BC0_type==0 .and. zeon(1,m) < psimin)then
                  zeon(1,m)=zeon0(1,m)
                  zeon(2,m)=2.0_wp*pi-zeon0(2,m)
                  ! Do not reset zeon(3,m) or toroidal angle !!
                  ! This is to ensure that this particle is in the same toroidal rank.
                  zeon(4,m)=zeon0(4,m)
                  zeon(5,m)=zeon0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zeon(5,m)=-zeon0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zeon(5,m)=0.0       ! ? more or less, close to delta-f=0 at inner boundary
                  zeon(7,m)=zeon0(7,m)
                  if(Flag_w_Evol>=3) zeon(8,m)=zeon0(8,m)
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


         if(Flag_Ptl_BC0_type==2 .or. Flag_Ptl_BC1_type==2) then
            ! Do not remove any particles here.
            ! (They will be removed by a separate routine every n-th timestep)
            ! Here, change the weight of original particles (labels>0) properly
            ! and let them flow around near the boundary.
            ! This method could be beneficial to fill the phase-space more efficiently.

            !$acc parallel loop present(zeon,zeon0,elabels,back_f0_e)
            do m=m0,m1
               ! If particles get into the region beyond boundaries
               ! reset their weights and labels properly while keeping their dxdv
               if(Flag_Ptl_BC1_type==2 .and. zeon(1,m)>psimax)then
                  ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
                  ! zeon0(9,m) is f_tot
                  ! TODO: check the definition of adiabterm outside the physical boundary
                  ! Here, assume that f_AD is smooth across the physical boundary so that zeon(9,m) is unchanged
                  integral_f0 = zeon0(6,m)*back_f0_e(m)
                  zeon(8,m) = integral_f0 -int(ifelectron/2)*zeon(9,m)! full-f (ifelectron<2) or full-h (ifelectron==2)
                  zeon(5,m) = zeon(8,m)-integral_f0 ! df or dh
                  zeon0(9,m) = integral_f0 ! f_tot

                  is_this_trc=0
                  !$acc loop seq
                  do k=1,num_of_etracers
                     if (elabels(m) == etracers_labels(k)) then
                        is_this_trc=1 ! this is one of tracers
                        exit
                     endif
                  enddo
                  if(is_this_trc==0) then
                     elabels(m) = -3.0 ! change elables value only if this is not a tracer
                  endif
               elseif(Flag_Ptl_BC0_type==2 .and. zeon(1,m)<psimin)then
                  ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
                  ! zeon0(9,m) is f_tot
                  ! TODO: check the definition of adiabterm outside the physical boundary
                  ! Here, assume that f_AD is smooth across the physical boundary so that zeon(9,m) is unchanged
                  integral_f0 = zeon0(6,m)*back_f0_e(m)
                  zeon(8,m) = integral_f0 -int(ifelectron/2)*zeon(9,m)! full-f (ifelectron<2) or full-h (ifelectron==2)
                  zeon(5,m) = zeon(8,m)-integral_f0 ! df or dh
                  zeon0(9,m) = integral_f0 ! f_tot

                  is_this_trc=0
                  !$acc loop seq
                  do k=1,num_of_etracers
                     if (elabels(m) == etracers_labels(k)) then
                        is_this_trc=1 ! this is one of tracers
                        exit
                     endif
                  enddo
                  if(is_this_trc==0) then
                     elabels(m) = -1.0 ! change elables value only if this is not a tracer
                  endif
               endif
            enddo
            acc_newer_D_zeon=.true.
         endif

      end subroutine treat_electrons_out_of_boundary


      subroutine treat_ions_out_of_boundary(m_start, m_end)

         use global_parameters
         use particle_array
         use field_array
         use function_gg_module
         use Ptl_locators_partial
         use Structured_Grid
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0, m1

         integer m,i,j,k,is_this_trc
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
         real(wp) target_w_tot, integral_f0

         real(wp) P_s,P_gq,P_gf,ws0,ws1,wgq0,wgq1,gyro_radius_in_s
         integer :: sid,gqid

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
                  zion(1,m)=zion0(1,m)
                  zion(2,m)=2.0_wp*pi-zion0(2,m)
                  ! Do not reset zion(3,m) or toroidal angle !!
                  ! This is to ensure that this particle is in the same toroidal rank.
                  zion(4,m)=zion0(4,m)
                  zion(5,m)=zion0(5,m) ! this sounds reasonable: energy, momentum and particle conserving boundary
                  !   zion(5,m)=-zion0(5,m)   ! remove density (weight) accumulation at edge when ifneoclassical=1
                  !   zion(5,m)=0.0       ! ? more or less, close to delta-f=0 at outer boundary
                  zion(7,m)=zion0(7,m)
               elseif(Flag_Ptl_BC0_type==0 .and. zion(1,m) < psimin)then
                  zion(1,m)=zion0(1,m)
                  zion(2,m)=2.0_wp*pi-zion0(2,m)
                  ! Do not reset zion(3,m) or toroidal angle !!
                  ! This is to ensure that this particle is in the same toroidal rank.
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

         if(Flag_Ptl_BC0_type==2 .or. Flag_Ptl_BC1_type==2) then
            ! Do not remove any particles here.
            ! (They will be removed by a separate routine every n-th timestep)
            ! Here, change the weight of original particles (labels>0) properly
            ! and let them flow around near the boundary.
            ! This method could be beneficial to fill the phase-space more efficiently.

            !$acc parallel loop present(zion,zion0,gyro_over_sqrt_mu,ilabels,back_f0_i)
            do m=m0,m1
               P_s = zion(7,m)
               P_gq = modulo(zion(2,m),2*pi)

               sid = max(1,min(int((P_s-s0)*SG%inv_ds) +1,SG%dim_s-1)) ! +1 is due to Fortran's default index starting from 1
               ws1 = (P_s-SG%s_1D(sid))*SG%inv_ds
               ws0 = (1.0-ws1)
               gqid = P_gq*SG%inv_dgq+1
               wgq1 = (P_gq-SG%gq_1D(gqid))*SG%inv_dgq
               wgq0 = 1.0-wgq1

               ! Calculate radial gyro radius in s (not a)
               gyro_radius_in_s  =  zion(6,m)*( ws0*(wgq0*gyro_over_sqrt_mu(1,1,gqid,sid)   + wgq1*gyro_over_sqrt_mu(1,1,gqid+1,sid)) &
                                               +ws1*(wgq0*gyro_over_sqrt_mu(1,1,gqid,sid+1) + wgq1*gyro_over_sqrt_mu(1,1,gqid+1,sid+1)))

               ! If particles get into the region beyond boundaries
               ! reset their weights and labels properly while keeping their dxdv
               if(Flag_Ptl_BC1_type==2) then
                  if(Ptl_BC1_val==0) then
                     ! If BDY1 is the wall condition, we should consider the gyroradius to determinte whether the gyro-orbit touches the wall
                     ! P_s = P_s+gyro_radius_in_s
                     if(zion(7,m)+gyro_radius_in_s>s1)then
                        integral_f0 = zion0(6,m)*back_f0_i(m) ! f_0
                        zion0(8,m) = 0.0 ! This is the wall condition, ftot=0 for the ion touching the wall (considering gyroradius)
                        zion(5,m) = -integral_f0 ! ftot=f0+df=0 --> df=-f0
                        zion(8,m) = integral_f0 ! f_0

                        is_this_trc=0
                        !$acc loop seq
                        do k=1,num_of_itracers
                           if (ilabels(m) == itracers_labels(k)) then
                              is_this_trc=1 ! this is one of tracers
                              exit
                           endif
                        enddo
                        if(is_this_trc==0) then
                        ilabels(m) = -3.0
                     endif
                     endif
                  else
                     if(zion(7,m)>s1)then
                        integral_f0 = zion0(6,m)*back_f0_i(m) ! f_0
                        zion0(8,m) = integral_f0 ! full-f
                        zion(5,m) = 0.0! df
                        zion(8,m) = integral_f0 ! f_0
                        
                        is_this_trc=0
                        !$acc loop seq
                        do k=1,num_of_itracers
                           if (ilabels(m) == itracers_labels(k)) then
                              is_this_trc=1 ! this is one of tracers
                              exit
                           endif
                        enddo
                        if(is_this_trc==0) then
                        ilabels(m) = -3.0
                     endif
                  endif
                  endif
               elseif(Flag_Ptl_BC0_type==2) then
                  if(Ptl_BC0_val==0) then
                     ! If BDY0 is the wall condition, we should consider the gyroradius to determinte whether the gyro-orbit touches the wall
                     if(zion(7,m)-gyro_radius_in_s<s0) then
                        integral_f0 = zion0(6,m)*back_f0_i(m) ! f_0
                        zion0(8,m) = 0.0 ! This is the wall condition, ftot=0 for the ion touching the wall (considering gyroradius)
                        zion(5,m) = -integral_f0 ! ftot=f0+df=0 --> df=-f0
                        zion(8,m) = integral_f0 ! f_0
                        
                        is_this_trc=0
                        !$acc loop seq
                        do k=1,num_of_itracers
                           if (ilabels(m) == itracers_labels(k)) then
                              is_this_trc=1 ! this is one of tracers
                              exit
                           endif
                        enddo
                        if(is_this_trc==0) then
                           ilabels(m) = -1.0
                        endif
                     endif
                  else
                     if(zion(7,m)<s0) then
                        integral_f0 = zion0(6,m)*back_f0_i(m) ! f_0
                        zion0(8,m) = integral_f0 ! full-f
                        zion(5,m) = 0.0! df
                        zion(8,m) = integral_f0 ! f_0

                        is_this_trc=0
                        !$acc loop seq
                        do k=1,num_of_itracers
                           if (ilabels(m) == itracers_labels(k)) then
                              is_this_trc=1 ! this is one of tracers
                              exit
                           endif
                        enddo
                        if(is_this_trc==0) then
                           ilabels(m) = -1.0
                        endif
                     endif
                  endif
               endif
            enddo
            acc_newer_D_zion=.true.
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
         integer :: num_of_recycling_ptls
         integer i,m,k,ke,ki
         integer me0,mi0
         real(wp) :: ErgDampRatio=1e-1
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
            zion0(6,m) = sqrt(ErgDampRatio)*zion0(6,m) ! since v_perp is reduced by factor of sqrt(ErgDampRatio), Jvperp=2*pi*vperp is also reduced
            zion0(8,m) = zion0(8,m)/(outION%Recycling_Ratio) ! divide by Recycling Ratio to statistically match 100% Recycling
            zion(5,m)  = zion0(8,m) - back_f0_i(m)*zion0(6,m)
            ! zion0(6,m) = 1e-6*zion0(6,m) ! Reduce the phase volume to reduce the particle noise
         enddo

         !$acc parallel loop present(zion,zion0,wi_init,wi_vpara,wi_vperp,wi_accel,wi_fdamp,UnNorm_theta_i,UnNorm_zeta_i)
         do m=mi0,mi
            wi_init(m)  = zion(5,m)
            wi_vpara(m) = 0
            wi_vperp(m) = 0
            wi_accel(m) = 0
            wi_fdamp(m) = 0

            !$acc loop seq
            do k=1,5
               zion0(k,m)     = zion(k,m)
            enddo
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

            ! Phase space volume
            ! The power 1.5 is because of (v_abs^2) in area of velocity domain and (v_perp) of velocity Jacobian
            zeon0(6,ke)=zion0(6,ki)*(backTe(ke)/backT(ki))**1.5
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

            !$acc loop seq
            do k=1,5
               zeon0(k,m)     = zeon(k,m)
            enddo
            zeon0(7,m)       = zeon(7,m)
            UnNorm_theta_e(m)= zeon(2,m)
            UnNorm_zeta_e(m) = zeon(3,m)
         enddo


         call Clear_OutParticles

         deallocate(rand_arr)

      end function Produce_Recycling_Ion_and_Electrons_using_OutParticles


      subroutine Generate_Ions_beyond_VelocityCutOff_Fixed_Grid
         use global_parameters
         use particle_decomp
         use Structured_Grid
         use mpi_defs
         implicit none

         real(wp), allocatable :: tmp_a(:),tmp_gq(:),tmp_gf(:),tmp_upara(:),tmp_uperp(:)
         real(wp), allocatable :: mesh_norm_dxdv(:,:,:,:,:)
         real(wp) den
         real(wp) upara_min, upara_max, uperp_min, uperp_max
         integer dim_upara, dim_uperp
         integer Nm,i, loop
         real(wp) time0,time1
         real(wp) upara_min_arr(3), upara_max_arr(3)
         real(wp) uperp_min_arr(3), uperp_max_arr(3)

         den = ini_mi/((a1-a0)*(2*pi)*(zetamax-zetamin)*(2*umax)*umax)

         upara_min_arr(1:3) = [-1.1*umax, -umax,     umax ]
         upara_max_arr(1:3) = [    -umax, +umax, 1.1*umax ]

         uperp_min_arr(1:3) = [0.0,            umax, 0.0    ]
         uperp_max_arr(1:3) = [ 1.1*umax, 1.1*umax, 1.1*umax]

         do loop=1,size(upara_min_arr)

            ! upara_min = 2.0; upara_max = 5.0; dim_upara = 5
            ! uperp_min = 1.0; uperp_max = 4.0; dim_uperp = 5
            upara_min = upara_min_arr(loop); upara_max = upara_max_arr(loop); dim_upara = 6
            uperp_min = uperp_min_arr(loop); uperp_max = uperp_max_arr(loop); dim_uperp = 6


            allocate(mesh_norm_dxdv(dim_upara,dim_uperp,PSCV%dim_a,PSCV%dim_gq,2))
            allocate(tmp_upara(dim_upara),tmp_uperp(dim_uperp),tmp_gf(2))

            tmp_gf(1) = zetamin
            tmp_gf(2) = zetamax
            do i=1,dim_upara
               tmp_upara(i) = upara_min + (i-1)*(upara_max-upara_min)/(dim_upara-1)
            enddo
            do i=1,dim_uperp
               tmp_uperp(i) = uperp_min + (i-1)*(uperp_max-uperp_min)/(dim_uperp-1)
            enddo


            Nm = den*(PSCV%a_max-PSCV%a_min)*(2.0*pi)*(tmp_gf(2)-tmp_gf(1))*(upara_max-upara_min)*(uperp_max-uperp_min)
            if(mype==0) then
               write(stdout,'(3(A,I0))') "gen_SRC_vel Nm_i=",Nm,"; mi=", mi,"; mimax=",mimax
               call flush(stdout)
            endif

            call Measure_Global_mesh_norm_dxdv_of_Ions(tmp_upara, tmp_uperp, mesh_norm_dxdv)
            call Generate_Ions_in_given_Phase_Space(tmp_upara,tmp_uperp,PSCV%a_1D,PSCV%gq_1D,tmp_gf,Nm,mesh_norm_dxdv,-2.0)


            deallocate(tmp_upara,tmp_uperp,tmp_gf)
            deallocate(mesh_norm_dxdv)

         enddo !loop

      end subroutine Generate_Ions_beyond_VelocityCutOff_Fixed_Grid


      subroutine Generate_Ions_beyond_PhysicalBoundary
         use global_parameters
         use particle_decomp
         use mpi_defs
         implicit none

         real(wp), allocatable :: tmp_a(:),tmp_gq(:),tmp_gf(:),tmp_upara(:),tmp_uperp(:)
         real(wp) den
         integer Nm,i, loop


         den = ini_mi/((a1-a0)*(2*pi)*(zetamax-zetamin)*(2*umax)*umax)
         allocate(tmp_upara(2),tmp_uperp(2),tmp_a(2),tmp_gq(2),tmp_gf(2))

         if(Flag_Ptl_BC0_type==1 .or. Flag_Ptl_BC0_type==2) then
            tmp_upara = [-1.1*umax, 1.1*umax]
            tmp_uperp = [0.0, 1.1*umax]
            tmp_a = [a0-0.05, a0]
            tmp_gq = [0.0, 2.0*pi]
            tmp_gf = [zetamin, zetamax]

            Nm = den*(tmp_upara(2)-tmp_upara(1))*(tmp_uperp(2)-tmp_uperp(1))*(tmp_a(2)-tmp_a(1))*(tmp_gq(2)-tmp_gq(1))*(tmp_gf(2)-tmp_gf(1))
            if(mype==0) then
               write(stdout,'(3(A,I0))') "gen_SRC inner radial Nm_i=",Nm,"; mi=", mi,"; mimax=",mimax
               call flush(stdout)
            endif

            call Generate_Ions_in_given_Phase_Space(tmp_upara,tmp_uperp,tmp_a,tmp_gq,tmp_gf,Nm,in_label_val=-1.0)
         endif

         if(Flag_Ptl_BC1_type==1 .or. Flag_Ptl_BC1_type==2) then
            tmp_upara = [-umax, umax]
            tmp_uperp = [0.0, umax]
            tmp_a = [a1, a1+0.05]
            tmp_gq = [0.0, 2.0*pi]
            tmp_gf = [zetamin, zetamax]

            Nm = den*(tmp_upara(2)-tmp_upara(1))*(tmp_uperp(2)-tmp_uperp(1))*(tmp_a(2)-tmp_a(1))*(tmp_gq(2)-tmp_gq(1))*(tmp_gf(2)-tmp_gf(1))
            if(mype==0) then
               write(stdout,'(3(A,I0))') "gen_SRC outer radial Nm_i=",Nm,"; mi=", mi,"; mimax=",mimax
               call flush(stdout)
            endif

            call Generate_Ions_in_given_Phase_Space(tmp_upara,tmp_uperp,tmp_a,tmp_gq,tmp_gf,Nm,in_label_val=-3.0)
         endif

         deallocate(tmp_upara,tmp_uperp,tmp_a,tmp_gq,tmp_gf)
      end subroutine Generate_Ions_beyond_PhysicalBoundary

      subroutine Generate_Ions_beyond_Boundaries_using_Adaptive_PSCV
         use global_parameters
         use particle_decomp
         use Structured_Grid
         use mpi_defs
         implicit none

         real(wp), allocatable :: tmp_gf(:)
         real(wp), allocatable :: mesh_norm_dxdv(:,:,:,:,:), mesh_MarkerDen(:,:,:,:,:)
         integer dim_ur,dim_uq
         integer i,j
         real(wp), allocatable :: ur_bound(:,:,:,:)
         real(wp) PSCV_mean_a

         if(Flag_Load_and_WeightDef==1 .or. Flag_Load_and_WeightDef==2) then
            dim_ur=5
            dim_uq=8
            allocate(tmp_gf(2))
            allocate(mesh_norm_dxdv(dim_ur,dim_uq,PSCV%dim_a,PSCV%dim_gq,2), mesh_MarkerDen(dim_ur,dim_uq,PSCV%dim_a,PSCV%dim_gq,2))

            tmp_gf(1) = zetamin
            tmp_gf(2) = zetamax

            allocate(ur_bound(PSCV%dim_a,PSCV%dim_gq,size(tmp_gf),2))

            PSCV_mean_a = 0.5*(PSCV%a_min+PSCV%a_max)

            if( (PSCV_mean_a-a0)*(PSCV_mean_a-a1)<0) then
               ur_bound(:,:,:,1) = umax ! Inside the wall
            else
               ur_bound(:,:,:,1) = 0.0 ! Outside the wall
            endif

            do i=1,PSCV%dim_a
               do j=1,PSCV%dim_gq
                  ur_bound(i,j,:,2) = umax + PSCV%Ion_dupara_per_nstep_SM(j,i,:)
                  ! ur_bound(i,j,:,2) = umax + PSCV%Ion_dupara_per_nstep(j,i,:)
               enddo
            enddo

            ! call Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions([10,10],mesh_norm_dxdv)
            ! call Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions([10,10],mesh_MarkerDen)

            call Measure_Global_mesh_norm_dxdv_of_Ions_PolarVelCoord_on_TgtUr(dim_ur,dim_uq,ur_bound, mesh_norm_dxdv, mesh_MarkerDen, InFlag_debug_output=.false.)
            call Generate_Ions_in_given_Phase_Space_PolarVelCoord_AdapTgtUr(PSCV%a_1D,PSCV%gq_1D,tmp_gf,ur_bound, mesh5D_norm_occp_dxdv=mesh_norm_dxdv,in_label_val=PSCV%ptl_label)
            ! ! Using MarkerDen is possible, but its accuracy could be low if total number of particles are not enough
            ! call Generate_Ions_in_given_Phase_Space_PolarVelCoord_AdapTgtUr(PtlSrcCV%a_1D,PtlSrcCV%gq_1D,tmp_gf,ur_bound,mesh_occp_MarkerDen=mesh_MarkerDen, mesh5D_norm_occp_dxdv=mesh_norm_dxdv,in_label_val=-2.0)

            deallocate(ur_bound)
            deallocate(tmp_gf)
            deallocate(mesh_norm_dxdv,mesh_MarkerDen)
         endif


      end subroutine Generate_Ions_beyond_Boundaries_using_Adaptive_PSCV


      subroutine Generate_Electrons_beyond_Boundaries_using_Adaptive_PSCV
         use global_parameters
         use particle_decomp
         use Structured_Grid
         use mpi_defs
         implicit none

         real(wp), allocatable :: tmp_gf(:)
         real(wp), allocatable :: mesh_norm_dxdv(:,:,:,:,:), mesh_MarkerDen(:,:,:,:,:)
         integer dim_ur, dim_uq
         integer i,j
         real(wp), allocatable :: ur_bound(:,:,:,:)
         real(wp) PSCV_mean_a

         if(Flag_Load_and_WeightDef==1 .or. Flag_Load_and_WeightDef==2) then
            dim_ur=5
            dim_uq=8
            allocate(tmp_gf(2))
            allocate(mesh_norm_dxdv(dim_ur,dim_uq,PSCV%dim_a,PSCV%dim_gq,2), mesh_MarkerDen(dim_ur,dim_uq,PSCV%dim_a,PSCV%dim_gq,2))

            tmp_gf(1) = zetamin
            tmp_gf(2) = zetamax

            allocate(ur_bound(PSCV%dim_a,PSCV%dim_gq,size(tmp_gf),2))

            PSCV_mean_a = 0.5*(PSCV%a_min+PSCV%a_max)
            if( (PSCV_mean_a-a0)*(PSCV_mean_a-a1)<0) then
               ur_bound(:,:,:,1) = umax ! Inside the wall
            else
               ur_bound(:,:,:,1) = 0.0 ! Outside the wall
            endif

            do i=1,PSCV%dim_a
               do j=1,PSCV%dim_gq
                  ur_bound(i,j,:,2) = umax + PSCV%Elec_dupara_per_nstep_SM(j,i,:)
                  ! ur_bound(i,j,:,2) = umax + PSCV%Elec_dupara_per_nstep(j,i,:)
               enddo
            enddo


            call Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord_on_TgtUr(dim_ur,dim_uq,ur_bound, mesh_norm_dxdv, mesh_MarkerDen, InFlag_debug_output=.false.)
            ! call Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord_on_TgtUr(dim_ur,dim_uq,ur_bound, mesh_norm_dxdv, mesh_MarkerDen, InFlag_debug_output=.true.)
            ! call Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord_on_TgtUr(dim_ur,dim_uq,ur_bound, mesh_norm_dxdv, mesh_MarkerDen, InFlag_debug_output=.true.)
            call Generate_Electrons_in_given_Phase_Space_PolarVelCoord_AdapTgtUr(PSCV%a_1D,PSCV%gq_1D,tmp_gf,ur_bound, mesh5D_norm_occp_dxdv=mesh_norm_dxdv,in_label_val=PSCV%ptl_label)
            ! call Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord_on_TgtUr(dim_ur,dim_uq,ur_bound, mesh_norm_dxdv, mesh_MarkerDen, InFlag_debug_output=.true.)
            ! ! Using MarkerDen is possible, but its accuracy could be low if total number of particles are not enough
            ! call Generate_Electrons_in_given_Phase_Space_PolarVelCoord_AdapTgtUr(PtlSrcCV%a_1D,PtlSrcCV%gq_1D,tmp_gf,ur_bound,mesh_occp_MarkerDen=mesh_MarkerDen, mesh5D_norm_occp_dxdv=mesh_norm_dxdv,in_label_val=-2.0)

            deallocate(ur_bound)
            deallocate(tmp_gf)
            deallocate(mesh_norm_dxdv,mesh_MarkerDen)
         endif


      end subroutine Generate_Electrons_beyond_Boundaries_using_Adaptive_PSCV

      subroutine Generate_Electrons_beyond_PhysicalBoundary
         use global_parameters
         use particle_decomp
         use mpi_defs
         implicit none

         real(wp), allocatable :: tmp_a(:),tmp_gq(:),tmp_gf(:),tmp_upara(:),tmp_uperp(:),tmp_ur(:),tmp_uq(:)
         real(wp) ini_Ne_marker_density
         integer Nm,i, loop


         ini_Ne_marker_density = ini_me/((a1-a0)*(2*pi)*(zetamax-zetamin)*(2*umax)*umax)

         allocate(tmp_upara(2),tmp_uperp(2),tmp_a(2),tmp_gq(2),tmp_gf(2),tmp_ur(2),tmp_uq(2))

         if(Flag_Ptl_BC0_type==1 .or. Flag_Ptl_BC0_type==2) then
            tmp_a = [a0-0.05, a0]
            tmp_gq = [0.0, 2.0*pi]
            tmp_gf = [zetamin, zetamax]

            if(Flag_Load_and_WeightDef==1 .or. Flag_Load_and_WeightDef==2) then
               tmp_ur = [0.0, 1.1*umax]
               tmp_uq = [0.0, pi]
               Nm = ini_Ne_marker_density*(0.5*pi*(tmp_ur(2)**2.0-tmp_ur(1)**2.0))*(tmp_a(2)-tmp_a(1))*(tmp_gq(2)-tmp_gq(1))*(tmp_gf(2)-tmp_gf(1))
               call Generate_Electrons_in_given_Phase_Space_PolarVelCoord(tmp_ur,tmp_uq,tmp_a,tmp_gq,tmp_gf,Nm, in_label_val=-1.0)
            endif

            if(mype==0) then
               write(stdout,'(3(A,I0))') "gen_SRC inner radial Nm_e=",Nm,"; me=", me,"; memax=",memax
               call flush(stdout)
            endif
         endif

         if(Flag_Ptl_BC1_type==1 .or. Flag_Ptl_BC1_type==2) then
            tmp_a = [a1, a1+0.05]
            tmp_gq = [0.0, 2.0*pi]
            tmp_gf = [zetamin, zetamax]

            if(Flag_Load_and_WeightDef==1 .or. Flag_Load_and_WeightDef==2) then
               tmp_ur = [0.0, 1.1*umax]
               tmp_uq = [0.0, pi]
               Nm = ini_Ne_marker_density*(0.5*pi*(tmp_ur(2)**2.0-tmp_ur(1)**2.0))*(tmp_a(2)-tmp_a(1))*(tmp_gq(2)-tmp_gq(1))*(tmp_gf(2)-tmp_gf(1))
               call Generate_Electrons_in_given_Phase_Space_PolarVelCoord(tmp_ur,tmp_uq,tmp_a,tmp_gq,tmp_gf,Nm, in_label_val=-3.0)
            endif

            if(mype==0) then
               write(stdout,'(3(A,I0))') "gen_SRC outer radial Nm_e=",Nm,"; me=", me,"; memax=",memax
               call flush(stdout)
            endif
         endif

         deallocate(tmp_upara,tmp_uperp,tmp_a,tmp_gq,tmp_gf,tmp_ur,tmp_uq)
      end subroutine Generate_Electrons_beyond_PhysicalBoundary

      subroutine Measure_Global_mesh_norm_dxdv_of_Ions(upara,uperp,mesh_norm_dxdv)
         use global_parameters
         use particle_array
         use particle_decomp
         use Structured_Grid
         use esi_interface
         use mpi_defs
         use function_gg_module
         use hdf5
         use h5lt
         implicit none
         real(wp), intent(in), dimension(:) :: upara,uperp
         real(wp), intent(inout) :: mesh_norm_dxdv(:,:,:,:,:) ! 5-D grid (upara,uperp,a,gq,gf) normalized dxdv by pre-occupying particles
         real(wp), allocatable :: tmp_mesh(:,:,:,:,:)

         real(wp), allocatable :: a_1D(:), gq(:), gf(:)
         real(wp) :: mesh_min_x(5),mesh_max_x(5), mesh_inv_dx(5)
         integer :: midx(5), mdims(5) ! mesh index and dimensions
         real(wp) :: ptl_data(5)
         integer m,i,j,k,length,ierror,istatus(MPI_STATUS_SIZE)
         integer aid,gqid,gfid,uparaid,uperpid
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max, upara_min, upara_max, uperp_min, uperp_max
         real(wp) :: da, dgq, dgf, dupara, duperp
         integer :: dim_a, dim_gq, dim_gf, dim_upara, dim_uperp
         logical state_in_mesh ! true=in, false=out of given mesh domain
         real(wp), external :: temperature
         real(wp) mesh_cell_vol ! 5-D mesh volume
         real(wp) time0, time1

         character(len=128) filename
         integer(HID_T) file_id, gid
         integer :: h5err
         integer(HSIZE_T) :: h5dims(5)

         ! Global mesh means that the mesh covers all the spatial range (a0-a1, 0-2pi, 0-2pi)
         ! by taking into account all processors and their communications.
         ! Also, the periodical conditions into theta and toroidal direcitons are considered


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         ! ! For simplicity, let's use a zero-th order shape function (closest near grid)
         if(acc_newer_H_zion) then
            !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_H_zion=.false.
            acc_newer_D_zion=.false.
         endif

         call n_T_u_group   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         call esigetbf(zion,mi,ziondim)   ! get B-->BB, I-->F  12/07
         !$acc parallel loop present(backu,F,BB)
         do m=1,mi
            backu(m)=backu(m)*F(m)/BB(m) ! =u_//=I*omega_t/B
         enddo

         ! allocate a,gq,gf and copy SG information into
         allocate(a_1D(PSCV%dim_a),gq(PSCV%dim_gq),gf(2))
         a_1D=PSCV%a_1D
         gq=PSCV%gq_1D
         gf(1)=zetamin
         gf(2)=zetamax

         upara_min = minval(upara); upara_max = maxval(upara); dim_upara = size(upara)
         uperp_min = minval(uperp); uperp_max = maxval(uperp); dim_uperp = size(uperp)
         a_min = minval(a_1D); a_max = maxval(a_1D); dim_a = size(a_1D)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         dupara=(upara_max-upara_min)/(dim_upara-1)
         duperp=(uperp_max-uperp_min)/(dim_uperp-1)
         da=(a_max-a_min)/(dim_a-1)
         dgq=(gq_max-gq_min)/(dim_gq-1)
         dgf=(gf_max-gf_min)/(dim_gf-1)

         mesh_min_x = [upara_min, uperp_min, a_min, gq_min, gf_min]
         mesh_max_x = [upara_max, uperp_max, a_max, gq_max, gf_max]
         mesh_inv_dx = [1.0/dupara, 1.0/duperp, 1.0/da, 1.0/dgq, 1.0/dgf]
         mdims = [dim_upara, dim_uperp, dim_a, dim_gq, dim_gf]

         ! Measure normalized phase volume from ions that reside in mesh domain
         mesh_norm_dxdv=0.0
         !$acc parallel loop private(ptl_data,midx)
         do m=1,mi
            ptl_data(1) = (zion(4,m)*BB(m) -backu(m))/sqrt(backT(m)) ! upara
            ptl_data(2) = zion(6,m)*sqrt(2.0*BB(m)/backT(m)) ! uperp
            ptl_data(3) = zion(1,m) ! a
            ptl_data(4) = zion(2,m) ! gq
            ptl_data(5) = zion(3,m) ! gf

            state_in_mesh = .true.
            !$acc loop seq
            do i=1,5
               state_in_mesh = state_in_mesh .and. ((ptl_data(i)-mesh_min_x(i))*(ptl_data(i)-mesh_max_x(i))<=0)
            enddo

            if(state_in_mesh) then
               !$acc loop seq
               do i=1,5
                  midx(i) = ceiling( (ptl_data(i)-mesh_min_x(i))*mesh_inv_dx(i) + 0.5 ) ! round operation to find nearst grid
               enddo
               !$acc atomic update
               mesh_norm_dxdv(midx(1),midx(2),midx(3),midx(4),midx(5)) = mesh_norm_dxdv(midx(1),midx(2),midx(3),midx(4),midx(5)) + zion0(6,m)
            endif
         enddo



         ! Communicate with other processors
         allocate(tmp_mesh(mdims(1),mdims(2),mdims(3),mdims(4),mdims(5)))
         if(npartdom>1) then
            tmp_mesh = mesh_norm_dxdv
            length = product(mdims)
            call MPI_ALLREDUCE(tmp_mesh, mesh_norm_dxdv, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
         endif

         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically although they have different gqid)
         mesh_norm_dxdv(:,:,:,1,:) = mesh_norm_dxdv(:,:,:,1,:) + mesh_norm_dxdv(:,:,:,dim_gq,:)
         mesh_norm_dxdv(:,:,:,dim_gq,:) = mesh_norm_dxdv(:,:,:,1,:)

         !toroidal communication
         length = product(mdims(1:4))
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1, tmp_mesh(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         mesh_norm_dxdv(:,:,:,:,dim_gf) = mesh_norm_dxdv(:,:,:,:,dim_gf) + tmp_mesh(:,:,:,:,dim_gf)
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1, mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

         ! Normalize the results
         do gqid=1,PSCV%dim_gq
            do aid=1,PSCV%dim_a
               vthc = sqrt(temperature(a_1D(aid),gq(gqid)))
               do uperpid=1,dim_uperp
                  ! mesh_cell_vol = (cell spatial volume)*(dvpara)*integral(2*pi*vperp)_(v0-dv/2)_(v0+dv/2)
                  if(uperpid==1) then
                     ! integral(vperp)_v0_(v0+0.5dv) = 0.5*(v0*dv+0.25*dv^2)
                     mesh_cell_vol = PSCV%Cell_dV(gqid,aid)*(dupara*vthc)*2.0*pi*0.5*(uperp(uperpid)*duperp+0.25*(duperp**2.0))*(vthc**2.0)
                  elseif(uperpid==dim_uperp) then
                      ! integral(vperp)_(v0-0.5dv)_(v0) = 0.5*(v0*dv-0.25*dv^2)
                     mesh_cell_vol = PSCV%Cell_dV(gqid,aid)*(dupara*vthc)*2.0*pi*0.5*(uperp(uperpid)*duperp-0.25*(duperp**2.0))*(vthc**2.0)
                  else
                     ! integral(vperp)_(v0-0.5dv)_(v0+0.5dv) = (v0*dv)
                     mesh_cell_vol = PSCV%Cell_dV(gqid,aid)*(dupara*vthc)*2.0*pi*(uperp(uperpid)*duperp)*(vthc**2.0)
                  endif

                  ! Normalize using Cell_dxdv
                  mesh_norm_dxdv(:,uperpid,aid,gqid,:) = mesh_norm_dxdv(:,uperpid,aid,gqid,:)/mesh_cell_vol
                  ! Double values at vpara edges since their cell volumes are half of the normal volume
                  mesh_norm_dxdv(1,uperpid,aid,gqid,:) = 2.0*mesh_norm_dxdv(1,uperpid,aid,gqid,:)
                  mesh_norm_dxdv(dim_upara,uperpid,aid,gqid,:) = 2.0*mesh_norm_dxdv(dim_upara,uperpid,aid,gqid,:)
               enddo
            enddo
         enddo


         ! !Write into hdf5 file
         ! if(mype==0) then
         !    write(filename, '("mesh_norm_Ion_dxdv_",i5.5,".h5")') istep/msnap
         !    call h5open_f(h5err)
         !    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)

         !    h5dims(1) = dim_upara
         !    call h5ltmake_dataset_f(file_id, "upara",1,h5dims,H5T_NATIVE_DOUBLE, upara, h5err)
         !    h5dims(1) = dim_uperp
         !    call h5ltmake_dataset_f(file_id, "uperp",1,h5dims,H5T_NATIVE_DOUBLE, uperp, h5err)
         !    h5dims(1) = dim_a
         !    call h5ltmake_dataset_f(file_id, "a",1,h5dims,H5T_NATIVE_DOUBLE, a, h5err)
         !    h5dims(1) = dim_gq
         !    call h5ltmake_dataset_f(file_id, "gq",1,h5dims,H5T_NATIVE_DOUBLE, gq, h5err)
         !    h5dims(1) = dim_gf
         !    call h5ltmake_dataset_f(file_id, "gf",1,h5dims,H5T_NATIVE_DOUBLE, gf, h5err)

         !    h5dims(1:5) = mdims(1:5)
         !    call h5ltmake_dataset_f(file_id, "mesh_norm_dxdv",5,h5dims,H5T_NATIVE_DOUBLE, mesh_norm_dxdv, h5err)

         !    call h5fclose_f(file_id, h5err)
         !    call h5close_f(h5err)
         ! endif

         deallocate(a_1D,gq,gf)
         deallocate(tmp_mesh)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Global_mesh_norm_dxdv_of_Ions takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Measure_Global_mesh_norm_dxdv_of_Ions

      subroutine Measure_Global_mesh_norm_dxdv_of_Electrons(upara,uperp,mesh_norm_dxdv)
         use global_parameters
         use particle_array
         use particle_decomp
         use Structured_Grid
         use esi_interface
         use mpi_defs
         use function_gg_module
         use hdf5
         use h5lt
         implicit none
         real(wp), intent(in), dimension(:) :: upara,uperp
         real(wp), intent(inout) :: mesh_norm_dxdv(:,:,:,:,:) ! 5-D grid (upara,uperp,a,gq,gf) normalized dxdv by pre-occupying particles
         real(wp), allocatable :: tmp_mesh(:,:,:,:,:)

         real(wp), allocatable :: a_1D(:), gq(:), gf(:)
         real(wp) :: mesh_min_x(5),mesh_max_x(5), mesh_inv_dx(5)
         integer :: midx(5), mdims(5) ! mesh index and dimensions
         real(wp) :: ptl_data(5)
         integer m,i,j,k,length,ierror,istatus(MPI_STATUS_SIZE)
         integer aid,gqid,gfid,uparaid,uperpid
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max, upara_min, upara_max, uperp_min, uperp_max
         real(wp) :: da, dgq, dgf, dupara, duperp
         integer :: dim_a, dim_gq, dim_gf, dim_upara, dim_uperp
         logical state_in_mesh ! true=in, false=out of given mesh domain
         real(wp), external :: temperature_e
         real(wp) mesh_cell_vol ! 5-D mesh volume
         real(wp) time0,time1

         character(len=128) filename, gname
         integer(HID_T) file_id, gid
         integer :: h5err
         integer(HSIZE_T) :: h5dims(5)

         ! Global mesh means that the mesh covers all the spatial range (a0-a1, 0-2pi, 0-2pi)
         ! by taking into account all processors and their communications.
         ! Also, the periodical conditions into theta and toroidal direcitons are considered

#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         ! ! For simplicity, let's use a zero-th order shape function (closest near grid)
         if(acc_newer_H_zeon) then
            !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif

         call n_T_u_group_e   ! get Te, note that back_f0_e cannot be calculated here because it requires zeon(4) and zeon(6)
         call esigetb(zeon,me,zeondim)   ! get B-->BB

         ! allocate a,gq,gf and copy SG information into
         allocate(a_1D(PSCV%dim_a),gq(PSCV%dim_gq),gf(2))
         a_1D=PSCV%a_1D
         gq=PSCV%gq_1D
         gf(1)=zetamin
         gf(2)=zetamax

         upara_min = minval(upara); upara_max = maxval(upara); dim_upara = size(upara)
         uperp_min = minval(uperp); uperp_max = maxval(uperp); dim_uperp = size(uperp)
         a_min = minval(a_1D); a_max = maxval(a_1D); dim_a = size(a_1D)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         dupara=(upara_max-upara_min)/(dim_upara-1)
         duperp=(uperp_max-uperp_min)/(dim_uperp-1)
         da=(a_max-a_min)/(dim_a-1)
         dgq=(gq_max-gq_min)/(dim_gq-1)
         dgf=(gf_max-gf_min)/(dim_gf-1)

         mesh_min_x = [upara_min, uperp_min, a_min, gq_min, gf_min]
         mesh_max_x = [upara_max, uperp_max, a_max, gq_max, gf_max]
         mesh_inv_dx = [1.0/dupara, 1.0/duperp, 1.0/da, 1.0/dgq, 1.0/dgf]
         mdims = [dim_upara, dim_uperp, dim_a, dim_gq, dim_gf]

         ! Measure normalized phase volume from electrons that reside in mesh domain
         mesh_norm_dxdv=0.0
         !$acc parallel loop private(ptl_data,midx)
         do m=1,me
            ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
            ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))
            ptl_data(1) = zeon(4,m)*csc0*BB(m)*sqrt(m0ms/(msm0*backTe(m)))! upara
            ptl_data(2) = zeon(6,m)*sqrt(2.0*BB(m))/sqrt(msm0*backTe(m))  ! uperp
            ptl_data(3) = zeon(1,m) ! a
            ptl_data(4) = zeon(2,m) ! gq
            ptl_data(5) = zeon(3,m) ! gf

            state_in_mesh = .true.
            !$acc loop seq
            do i=1,5
               state_in_mesh = state_in_mesh .and. ((ptl_data(i)-mesh_min_x(i))*(ptl_data(i)-mesh_max_x(i))<=0)
            enddo

            if(state_in_mesh) then
               !$acc loop seq
               do i=1,5
                  midx(i) = ceiling( (ptl_data(i)-mesh_min_x(i))*mesh_inv_dx(i) + 0.5 ) ! round operation to find nearst grid
               enddo
               !$acc atomic update
               mesh_norm_dxdv(midx(1),midx(2),midx(3),midx(4),midx(5)) = mesh_norm_dxdv(midx(1),midx(2),midx(3),midx(4),midx(5)) + zeon0(6,m)
            endif
         enddo



         ! Communicate with other processors
         allocate(tmp_mesh(mdims(1),mdims(2),mdims(3),mdims(4),mdims(5)))
         if(npartdom>1) then
            tmp_mesh = mesh_norm_dxdv
            length = product(mdims)
            call MPI_ALLREDUCE(tmp_mesh, mesh_norm_dxdv, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
         endif

         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically although they have different gqid)
         mesh_norm_dxdv(:,:,:,1,:) = mesh_norm_dxdv(:,:,:,1,:) + mesh_norm_dxdv(:,:,:,dim_gq,:)
         mesh_norm_dxdv(:,:,:,dim_gq,:) = mesh_norm_dxdv(:,:,:,1,:)

         !toroidal communication
         length = product(mdims(1:4))
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1, tmp_mesh(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         mesh_norm_dxdv(:,:,:,:,dim_gf) = mesh_norm_dxdv(:,:,:,:,dim_gf) + tmp_mesh(:,:,:,:,dim_gf)
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1, mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

         ! Normalize the results
         do gqid=1,PSCV%dim_gq
            do aid=1,PSCV%dim_a
               vthc = sqrt(temperature_e(a_1D(aid),gq(gqid))) ! Note that temperature_e gives the electron velocity not the energy. so not necessary to consider mass here
               do uperpid=1,dim_uperp
                  ! mesh_cell_vol = (cell spatial volume)*(dvpara)*integral(2*pi*vperp)_(v0-dv/2)_(v0+dv/2)
                  if(uperpid==1) then
                     ! integral(vperp)_v0_(v0+0.5dv) = 0.5*(v0*dv+0.25*dv^2)
                     mesh_cell_vol = PSCV%Cell_dV(gqid,aid)*(dupara*vthc)*2.0*pi*0.5*(uperp(uperpid)*duperp+0.25*(duperp**2.0))*(vthc**2.0)
                  elseif(uperpid==dim_uperp) then
                      ! integral(vperp)_(v0-0.5dv)_(v0) = 0.5*(v0*dv-0.25*dv^2)
                     mesh_cell_vol = PSCV%Cell_dV(gqid,aid)*(dupara*vthc)*2.0*pi*0.5*(uperp(uperpid)*duperp-0.25*(duperp**2.0))*(vthc**2.0)
                  else
                     ! integral(vperp)_(v0-0.5dv)_(v0+0.5dv) = (v0*dv)
                     mesh_cell_vol = PSCV%Cell_dV(gqid,aid)*(dupara*vthc)*2.0*pi*(uperp(uperpid)*duperp)*(vthc**2.0)
                  endif

                  ! Normalize using Cell_dxdv
                  mesh_norm_dxdv(:,uperpid,aid,gqid,:) = mesh_norm_dxdv(:,uperpid,aid,gqid,:)/mesh_cell_vol
                  ! Double values at vpara edges since their cell volumes are half of the normal volume
                  mesh_norm_dxdv(1,uperpid,aid,gqid,:) = 2.0*mesh_norm_dxdv(1,uperpid,aid,gqid,:)
                  mesh_norm_dxdv(dim_upara,uperpid,aid,gqid,:) = 2.0*mesh_norm_dxdv(dim_upara,uperpid,aid,gqid,:)

                  ! tmp_mesh(:,uperpid,aid,gqid,:) = mesh_cell_vol
                  ! tmp_mesh(1,uperpid,aid,gqid,:) = 0.5*mesh_cell_vol
                  ! tmp_mesh(dim_upara,uperpid,aid,gqid,:) = 0.5*mesh_cell_vol

               enddo
            enddo
         enddo


         ! !Write into hdf5 file
         ! if(mype==0) then
         !    ! write(filename, '("mesh_norm_Elec_dxdv_",i5.5,"_uparamin=",F0.2,".h5")') istep/msnap, upara_min
         !    write(filename, '("mesh_norm_Elec_dxdv_",i5.5,".h5")') istep/msnap

         !    call h5open_f(h5err)
         !    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err)
         !    if(h5err==-1) then ! If there is no file
         !       call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)
         !    endif

         !    if(0.5*(upara_min+upara_max)<-umax) then
         !       i=1
         !    elseif(0.5*(upara_min+upara_max)>umax) then
         !       i=2
         !    else
         !       i=3
         !    endif
         !    write(gname,'("case",I0)') i
         !    call h5gcreate_f(file_id, gname, gid, h5err)

         !    h5dims(1) = dim_upara
         !    call h5ltmake_dataset_f(gid, "upara",1,h5dims,H5T_NATIVE_DOUBLE, upara, h5err)
         !    h5dims(1) = dim_uperp
         !    call h5ltmake_dataset_f(gid, "uperp",1,h5dims,H5T_NATIVE_DOUBLE, uperp, h5err)
         !    h5dims(1) = dim_a
         !    call h5ltmake_dataset_f(gid, "a",1,h5dims,H5T_NATIVE_DOUBLE, a, h5err)
         !    h5dims(1) = dim_gq
         !    call h5ltmake_dataset_f(gid, "gq",1,h5dims,H5T_NATIVE_DOUBLE, gq, h5err)
         !    h5dims(1) = dim_gf
         !    call h5ltmake_dataset_f(gid, "gf",1,h5dims,H5T_NATIVE_DOUBLE, gf, h5err)

         !    h5dims(1:5) = mdims(1:5)
         !    call h5ltmake_dataset_f(gid, "mesh_norm_dxdv",5,h5dims,H5T_NATIVE_DOUBLE, mesh_norm_dxdv, h5err)
         !    call h5ltmake_dataset_f(gid, "mesh_cell_vol",5,h5dims,H5T_NATIVE_DOUBLE, tmp_mesh, h5err)

         !    call h5gclose_f(gid, h5err)
         !    call h5fclose_f(file_id, h5err)
         !    call h5close_f(h5err)
         ! endif

         deallocate(a_1D,gq,gf)
         deallocate(tmp_mesh)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Global_mesh_norm_dxdv_of_Electrons takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Measure_Global_mesh_norm_dxdv_of_Electrons


      subroutine Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord(ur,uq,mesh_norm_dxdv,InFlag_debug_output)
         use global_parameters
         use particle_array
         use particle_decomp
         use Structured_Grid
         use esi_interface
         use mpi_defs
         use function_gg_module
         use hdf5
         use h5lt
         implicit none
         real(wp), intent(in), dimension(:) :: ur,uq ! ur(radial), uq(theta)
         real(wp), intent(inout) :: mesh_norm_dxdv(:,:,:,:,:) ! 5-D grid (ur,uq,a,gq,gf) normalized dxdv by pre-occupying particles
         logical, intent(in), optional :: InFlag_debug_output
         logical :: flag_debug_output
         !$acc declare create(mesh_norm_dxdv)
         real(wp), allocatable :: tmp_mesh(:,:,:,:,:)
         ! real(wp), allocatable :: cell_vol_5D(:,:,:,:,:)

         real(wp) upara,uperp, Cell_dVol_x3, Cell_dVol_ur, Cell_dVol_uq
         real(wp), allocatable :: upara_2D(:,:), uperp_2D(:,:)
         real(wp), allocatable :: a_1D(:), gq(:), gf(:)
         real(wp) :: mesh_min_x(5),mesh_max_x(5), mesh_inv_dx(5)
         integer :: midx(5), mdims(5) ! mesh index and dimensions
         real(wp) :: ptl_data(5)
         integer m,i,j,k,LL,length,ierror,istatus(MPI_STATUS_SIZE)
         integer aid,gqid,gfid,urid,uqid
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max, ur_min, ur_max, uq_min, uq_max
         real(wp) :: da, dgq, dgf, dur, duq
         integer :: dim_a, dim_gq, dim_gf, dim_ur, dim_uq
         logical state_in_mesh ! true=in, false=out of given mesh domain
         real(wp), external :: temperature_e
         real(wp) mesh_cell_vol ! 5-D mesh volume
         real(wp) time0,time1
         real(wp), allocatable :: Te0_1D(:), Cell_dV(:,:)

         character(len=128) filename, gname
         integer(HID_T) file_id, gid
         integer :: h5err
         integer(HSIZE_T) :: h5dims(5)

         ! Global mesh means that the mesh covers all the spatial range (a0-a1, 0-2pi, 0-2pi)
         ! by taking into account all processors and their communications.
         ! Also, the periodical conditions into theta and toroidal direcitons are considered

#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         flag_debug_output = .false.
         if(present(InFlag_debug_output)) flag_debug_output=InFlag_debug_output

         ! ! For simplicity, let's use a zero-th order shape function (closest near grid)
         if(acc_newer_H_zeon) then
            !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif

         call n_T_u_group_e   ! get Te, note that back_f0_e cannot be calculated here because it requires zeon(4) and zeon(6)
         call esigetb(zeon,me,zeondim)   ! get B-->BB



         ! allocate a,gq,gf and copy SG information into
         allocate(a_1D(PSCV%dim_a),gq(PSCV%dim_gq),gf(2))
         a_1D=PSCV%a_1D
         gq=PSCV%gq_1D
         gf(1)=zetamin
         gf(2)=zetamax

         ur_min = minval(ur); ur_max = maxval(ur); dim_ur = size(ur)
         uq_min = minval(uq); uq_max = maxval(uq); dim_uq = size(uq)
         a_min = minval(a_1D); a_max = maxval(a_1D); dim_a = size(a_1D)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         dur=(ur_max-ur_min)/(dim_ur-1)
         duq=(uq_max-uq_min)/(dim_uq-1)
         da=(a_max-a_min)/(dim_a-1)
         dgq=(gq_max-gq_min)/(dim_gq-1)
         dgf=(gf_max-gf_min)/(dim_gf-1)

         mesh_min_x = [ur_min, uq_min, a_min, gq_min, gf_min]
         mesh_max_x = [ur_max, uq_max, a_max, gq_max, gf_max]
         mesh_inv_dx = [1.0/dur, 1.0/duq, 1.0/da, 1.0/dgq, 1.0/dgf]
         mdims = [dim_ur, dim_uq, dim_a, dim_gq, dim_gf]


         ! Measure normalized phase volume from electrons that reside in mesh domain
         ! mesh_norm_dxdv was allocated at the outside of this subroutine
         mesh_norm_dxdv=0.0
         !$acc update device(mesh_norm_dxdv)

         !$acc parallel loop present(mesh_norm_dxdv) private(ptl_data,midx)
         do m=1,me
            ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
            ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))
            upara = zeon(4,m)*csc0*BB(m)*sqrt(m0ms/(msm0*backTe(m)))! upara
            uperp = zeon(6,m)*sqrt(2.0*BB(m))/sqrt(msm0*backTe(m))  ! uperp

            ptl_data(1) = sqrt(upara**2.0+uperp**2.0) ! ur (raidus in the polar coordinate)
            ptl_data(2) = acos(upara/ptl_data(1)) ! uq (angle in the polar coordinate)
            ptl_data(3) = zeon(1,m) ! a
            ptl_data(4) = zeon(2,m) ! gq
            ptl_data(5) = zeon(3,m) ! gf

            state_in_mesh = .true.
            !$acc loop seq
            do i=1,5
               state_in_mesh = state_in_mesh .and. ((ptl_data(i)-mesh_min_x(i))*(ptl_data(i)-mesh_max_x(i))<=0)
            enddo

            if(state_in_mesh) then
               !$acc loop seq
               do i=1,5
                  midx(i) = ceiling( (ptl_data(i)-mesh_min_x(i))*mesh_inv_dx(i) + 0.5 ) ! round operation to find nearst grid
               enddo
               !$acc atomic update
               mesh_norm_dxdv(midx(1),midx(2),midx(3),midx(4),midx(5)) = mesh_norm_dxdv(midx(1),midx(2),midx(3),midx(4),midx(5)) + zeon0(6,m)
            endif
         enddo


         allocate(Te0_1D(PSCV%dim_a),Cell_dV(PSCV%dim_gq,PSCV%dim_a))

         Te0_1D = PSCV%Te0_1D
         Cell_dV = PSCV%Cell_dV

         ! allocate(cell_vol_5D(mdims(1),mdims(2),mdims(3),mdims(4),mdims(5)))

         length=dim_gq*dim_a*dim_uq*dim_ur
         ! Normalize the results
         !$acc parallel loop gang vector present(mesh_norm_dxdv) private(LL,k)
         do k=1,length
            LL=k
            gqid  = int((LL-1)/(dim_ur*dim_uq*dim_a))+1
            LL    = modulo((LL-1),dim_ur*dim_uq*dim_a)+1
            aid   = int((LL-1)/(dim_ur*dim_uq))+1
            LL    = modulo((LL-1),dim_ur*dim_uq)+1
            uqid  = int((LL-1)/(dim_ur))+1
            LL    = modulo((LL-1),dim_ur)+1
            urid  = LL

            vthc = sqrt(Te0_1D(aid)*m0ms)
               ! mesh_cell_vol = (cell spatial volume)*integral(2*pi*ur^2*sin(uq))
               !               = (cell spatial volume)*2*pi*integral(ur^2)_(ur0-dur/2)^(ur0+dur/2) * integral(sin(uq))_(uq0-duq/2)^(uq0+duq/2)
            ! Cell_dVol_x3 = PtlSrcCV%Cell_dV(gqid,aid)
            Cell_dVol_x3 = Cell_dV(gqid,aid)

                  if(uqid==1) then
                     ! integral(sin(uq))_(uq0)^(uq0+duq/2) = 2*sin(duq/4)*sin(uq0+duq/4)
                     Cell_dVol_uq = 2.0*sin(0.25*duq)*sin(uq(uqid)+0.25*duq)
                  elseif(uqid==dim_uq) then
                     ! integral(sin(uq))_(uq0-duq/2)^(uq0) = 2*sin(duq/4)*sin(uq0-duq/4)
                     Cell_dVol_uq = 2.0*sin(0.25*duq)*sin(uq(uqid)-0.25*duq)
                  else
                     ! integral(sin(uq))_(uq0-duq/2)^(uq0+duq/2) = 2*sin(duq/2)*sin(uq0)
                     Cell_dVol_uq = 2.0*sin(0.5*duq)*sin(uq(uqid))
                  endif

                     if(urid==1) then
                        ! integral(ur^2)_(ur0)_(ur0+0.5dur) = ((ur0+0.5dur)^3-(ur0)^3)/3
                        Cell_dVol_ur = ((ur(urid)+0.5*dur)**3.0-ur(urid)**3.0)/3.0
                     elseif(urid==dim_ur) then
                        ! integral(ur^2)_(ur0-0.5dur)_(ur0) = (ur0^3-(ur0-0.5dur)^3)/3
                        Cell_dVol_ur = (ur(urid)**3.0-(ur(urid)-0.5*dur)**3.0)/3.0
                     else
                        ! integral(ur^2)_(ur0-0.5dur)_(ur0+0.5dur) = ((ur0+0.5dur)^3-(ur0-0.5dur)^3)/3
                        Cell_dVol_ur = ((ur(urid)+0.5*dur)**3.0-(ur(urid)-0.5*dur)**3.0)/3.0
                     endif

                     mesh_cell_vol = 2.0*pi*(vthc**3.0)*Cell_dVol_x3*Cell_dVol_ur*Cell_dVol_uq

                     ! Normalize using Cell_dxdv
            !$acc loop seq
            do gfid=1,2
               mesh_norm_dxdv(urid,uqid,aid,gqid,gfid) = mesh_norm_dxdv(urid,uqid,aid,gqid,gfid)/mesh_cell_vol
               ! cell_vol_5D(urid,uqid,aid,gqid,gfid) = mesh_cell_vol
            enddo
         enddo
         deallocate(Te0_1D,Cell_dV)

         !$acc update self(mesh_norm_dxdv)

         allocate(tmp_mesh(mdims(1),mdims(2),mdims(3),mdims(4),mdims(5)))
         ! Communicate with other processors
         if(npartdom>1) then
            tmp_mesh = mesh_norm_dxdv
            length = product(mdims)
            call MPI_ALLREDUCE(tmp_mesh, mesh_norm_dxdv, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
         endif

         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically although they have different gqid)
         mesh_norm_dxdv(:,:,:,1,:) = mesh_norm_dxdv(:,:,:,1,:) + mesh_norm_dxdv(:,:,:,dim_gq,:)
         mesh_norm_dxdv(:,:,:,dim_gq,:) = mesh_norm_dxdv(:,:,:,1,:)

         !toroidal communication
         length = product(mdims(1:4))
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1, tmp_mesh(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         mesh_norm_dxdv(:,:,:,:,dim_gf) = mesh_norm_dxdv(:,:,:,:,dim_gf) + tmp_mesh(:,:,:,:,dim_gf)
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1, mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)



         !Write into hdf5 file
         if(mype==0 .and. flag_debug_output==.true.) then
            ! write(filename, '("polar_mesh_norm_Elec_dxdv_",i5.5,"_urmin=",F0.2,".h5")') istep/msnap, ur_min
            write(filename, '("mesh_norm_Elec_dxdv_",i5.5,".h5")') istep/msnap

            call h5open_f(h5err)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err)
            if(h5err==-1) then ! If there is no file
               call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)
            else ! If there is file
               write(filename, '("mesh_norm_Elec_after",i5.5,".h5")') istep/msnap
               call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)
            endif

            ! write(gname,'("polar")')
            ! call h5gcreate_f(file_id, gname, gid, h5err)
            gid=file_id

            h5dims(1) = dim_ur
            call h5ltmake_dataset_f(gid, "ur",1,h5dims,H5T_NATIVE_DOUBLE, ur, h5err)
            h5dims(1) = dim_uq
            call h5ltmake_dataset_f(gid, "uq",1,h5dims,H5T_NATIVE_DOUBLE, uq, h5err)
            h5dims(1) = dim_a
            call h5ltmake_dataset_f(gid, "a",1,h5dims,H5T_NATIVE_DOUBLE, a_1D, h5err)
            h5dims(1) = dim_gq
            call h5ltmake_dataset_f(gid, "gq",1,h5dims,H5T_NATIVE_DOUBLE, gq, h5err)
            h5dims(1) = dim_gf
            call h5ltmake_dataset_f(gid, "gf",1,h5dims,H5T_NATIVE_DOUBLE, gf, h5err)

            allocate(upara_2D(dim_ur,dim_uq),uperp_2D(dim_ur,dim_uq))
            do urid=1,dim_ur
               do uqid=1,dim_uq
                  upara_2D(urid,uqid) = ur(urid)*cos(uq(uqid))
                  uperp_2D(urid,uqid) = ur(urid)*sin(uq(uqid))
               enddo
            enddo

            h5dims(1:2) = [dim_ur,dim_uq]
            call h5ltmake_dataset_f(gid, "upara_2D",2,h5dims,H5T_NATIVE_DOUBLE, upara_2D, h5err)
            call h5ltmake_dataset_f(gid, "uperp_2D",2,h5dims,H5T_NATIVE_DOUBLE, uperp_2D, h5err)

            h5dims(1:5) = mdims(1:5)
            call h5ltmake_dataset_f(gid, "mesh_norm_dxdv",5,h5dims,H5T_NATIVE_DOUBLE, mesh_norm_dxdv, h5err)
            ! call h5ltmake_dataset_f(gid, "mesh_cell_vol",5,h5dims,H5T_NATIVE_DOUBLE, cell_vol_5D, h5err)

            tmp_mesh = mesh_norm_dxdv
            call Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions([10,10],tmp_mesh)
            call h5ltmake_dataset_f(gid, "mesh_norm_dxdv_SM",5,h5dims,H5T_NATIVE_DOUBLE, tmp_mesh, h5err)


            call h5gclose_f(gid, h5err)
            call h5fclose_f(file_id, h5err)
            call h5close_f(h5err)

            deallocate(upara_2D,uperp_2D)
         endif

         deallocate(a_1D,gq,gf)
         deallocate(tmp_mesh)
         ! deallocate(cell_vol_5D)
#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Global_mesh_norm_dxdv_of_Electrons_PolarCoord takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord

      subroutine Measure_Global_mesh_norm_dxdv_of_Ions_PolarVelCoord_on_TgtUr(dim_sr, dim_sq, tgt_ur_bound, mesh_norm_dxdv, mesh_MarkerDen, InFlag_debug_output)
         use global_parameters
         use particle_array
         use particle_decomp
         use Structured_Grid
         use esi_interface
         use mpi_defs
         use function_gg_module
         use hdf5
         use h5lt
         implicit none
         integer, intent(in) :: dim_sr, dim_sq
         real(wp), intent(in) :: tgt_ur_bound(:,:,:,:) ! (dim_a,dim_gq,dim_gf,2); upper ur and lower ur bound on 3-D grid (a,gq,gf)
         real(wp), intent(inout) :: mesh_norm_dxdv(:,:,:,:,:) ! 5-D grid (sr,sq,a,gq,gf) normalized dxdv by pre-occupying particles
         real(wp), intent(inout) :: mesh_MarkerDen(:,:,:,:,:) ! 5-D grid (sr,sq,a,gq,gf) normalized density by pre-occupying particles (on upace, no Jacobian)
         logical, intent(in), optional :: InFlag_debug_output
         logical :: flag_debug_output
         !$acc declare create(mesh_norm_dxdv,mesh_MarkerDen)
         real(wp), allocatable :: tmp_mesh(:,:,:,:,:),cell_vol_5D(:,:,:,:,:), cell_area_5D(:,:,:,:,:)
         real(wp), allocatable :: sr(:), sq(:)
         real(wp), allocatable :: a_1D(:), gq(:), gf(:)
         integer :: tmp_aid, tgt_aid, tgt_gqid, tgt_gfid, tgt_srid, tgt_sqid
         real(wp) ptl_ur, ptl_sr, ptl_sq, ptl_ur_bound(2)
         real(wp) tgt_gq, ws0,ws1, P_safety_q
         real(wp) upara,uperp, Cell_dVol_x3, Cell_dVol_sr, Cell_dVol_sq
         real(wp), allocatable :: spara_2D(:,:), sperp_2D(:,:), ur_4D(:,:,:,:), upara_5D(:,:,:,:,:), uperp_5D(:,:,:,:,:)
         integer m,i,j,k,LL,length,ierror,istatus(MPI_STATUS_SIZE)
         integer aid,gqid,gfid,srid,sqid
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max, sr_min, sr_max, sq_min, sq_max
         real(wp) :: da, dgq, dgf, dsr, dsq, inv_da, inv_dgq, inv_dgf, inv_dsr, inv_dsq
         real(wp) :: sr0,sr1,sq0,sq1,ur0,ur1
         integer :: dim_a, dim_gq, dim_gf
         real(wp) mesh_cell_vol, mesh_cell_area ! 5-D mesh volume and area of each mesh grid
         real(wp) time0,time1
         real(wp), allocatable :: Ti0_1D(:), Cell_dV(:,:)

         character(len=128) filename, gname
         integer(HID_T) file_id, gid
         integer :: h5err
         integer(HSIZE_T) :: h5dims(5)

         ! Global mesh means that the mesh covers all the spatial range (a0-a1, 0-2pi, 0-2pi)
         ! by taking into account all processors and their communications.
         ! Also, the periodical conditions into theta and toroidal direcitons are considered

#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         flag_debug_output = .false.
         if(present(InFlag_debug_output)) flag_debug_output=InFlag_debug_output

         ! ! For simplicity, let's use a zero-th order shape function (closest near grid)
         if(acc_newer_H_zion) then
            !$acc update device(zion(:,:mi),zion0(:,:mi),RK_cum_zion(:,:mi))
            acc_newer_H_zion=.false.
            acc_newer_D_zion=.false.
         endif

         call n_T_u_group   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         call esigetb(zion,mi,ziondim)   ! get B-->BB



         ! allocate a,gq,gf and copy SG information into
         allocate(a_1D(PSCV%dim_a),gq(PSCV%dim_gq),gf(2))
         allocate(sr(dim_sr),sq(dim_sq))
         a_1D=PSCV%a_1D
         gq=PSCV%gq_1D
         gf(1)=zetamin
         gf(2)=zetamax

         sr_min = 0; sr_max = 1;
         sq_min = 0; sq_max = pi;
         dsr=(sr_max-sr_min)/(dim_sr-1); inv_dsr = 1.0/dsr
         dsq=(sq_max-sq_min)/(dim_sq-1); inv_dsq = 1.0/dsq
         do i=1,dim_sr
            sr(i) = sr_min + (i-1)*dsr
         enddo
         do i=1,dim_sq
            sq(i) = sq_min + (i-1)*dsq
         enddo


         a_min = minval(a_1D); a_max = maxval(a_1D); dim_a = size(a_1D)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         da=(a_max-a_min)/(dim_a-1);     inv_da  = 1.0/da
         dgq=(gq_max-gq_min)/(dim_gq-1); inv_dgq = 1.0/dgq
         dgf=(gf_max-gf_min)/(dim_gf-1); inv_dgf = 1.0/dgf

         ! Measure normalized phase volume from electrons that reside in mesh domain
         ! mesh_norm_dxdv was allocated at the outside of this subroutine
         mesh_norm_dxdv=0.0
         mesh_MarkerDen=0.0
         !$acc update device(mesh_norm_dxdv,mesh_MarkerDen)

         !$acc parallel loop present(mesh_norm_dxdv,mesh_MarkerDen) private(ptl_ur_bound)
         do m=1,mi

            if( (zion(1,m)-a_min)*(zion(1,m)-a_max)>0 ) then
               ! If the particle is outside the PSCV domain, cycle the loop
               cycle
            endif

            tgt_aid  = ceiling((zion(1,m)-a_min)*inv_da+0.5) ! Find nearst grid
            tgt_gfid = ceiling((zion(3,m)-gf_min)*inv_dgf+0.5) ! Find nearst grid

            ! Interpolate safety factor
            tmp_aid = max(1,min(int((zion(1,m)-PSCV%a_min)*PSCV%inv_da) +1,PSCV%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
            ws1 = (zion(1,m)-PSCV%a_1D(tmp_aid))*PSCV%inv_da
            ws0 = (1.0-ws1)
            P_safety_q = ws0*PSCV%safety_q_1D(tmp_aid) + ws1*PSCV%safety_q_1D(tmp_aid+1)

            tgt_gq = modulo(zion(2,m) + (gf(tgt_gfid)-zion(3,m))/P_safety_q, 2.0*pi)
            tgt_gqid = ceiling((tgt_gq-gq_min)*inv_dgq+0.5) ! Find nearst grid

            ! Find ur bound of the particle
            ptl_ur_bound(1) = tgt_ur_bound(tgt_aid,tgt_gqid,tgt_gfid,1)
            ptl_ur_bound(2) = tgt_ur_bound(tgt_aid,tgt_gqid,tgt_gfid,2)

            ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
            ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))
            upara = zion(4,m)*BB(m)/sqrt(backT(m)) ! upara
            uperp = zion(6,m)*sqrt(2.0*BB(m))/sqrt(backT(m)) ! uperp

            ptl_ur = sqrt(upara**2.0+uperp**2.0)

            if( (ptl_ur-ptl_ur_bound(1))*(ptl_ur-ptl_ur_bound(2)) <0) then
               ! Measure particle's contribution to the norm_phase_volume and MarkerDen
               ptl_sr = (ptl_ur-ptl_ur_bound(1))/(ptl_ur_bound(2)-ptl_ur_bound(1))
               tgt_srid = ceiling((ptl_sr-sr_min)*inv_dsr+0.5)
               ptl_sq = acos(upara/ptl_ur)
               tgt_sqid = ceiling((ptl_sq-sq_min)*inv_dsq+0.5)

               !$acc atomic update
               mesh_norm_dxdv(tgt_srid,tgt_sqid,tgt_aid,tgt_gqid,tgt_gfid) = mesh_norm_dxdv(tgt_srid,tgt_sqid,tgt_aid,tgt_gqid,tgt_gfid) + zion0(6,m)
               !$acc atomic update
               mesh_MarkerDen(tgt_srid,tgt_sqid,tgt_aid,tgt_gqid,tgt_gfid) = mesh_MarkerDen(tgt_srid,tgt_sqid,tgt_aid,tgt_gqid,tgt_gfid) + 1.0
            endif
         enddo


         allocate(Ti0_1D(PSCV%dim_a),Cell_dV(PSCV%dim_gq,PSCV%dim_a))

         Ti0_1D = PSCV%Ti0_1D
         Cell_dV = PSCV%Cell_dV



         allocate(cell_vol_5D(dim_sr,dim_sq,dim_a,dim_gq,dim_gf))
         allocate(cell_area_5D(dim_sr,dim_sq,dim_a,dim_gq,dim_gf))

         length=dim_sq*dim_sr*dim_a*dim_gq*dim_gf
         ! Normalize the results
         !$acc parallel loop gang vector present(mesh_norm_dxdv,mesh_MarkerDen) private(LL,k,ptl_ur_bound)
         do k=1,length
            LL=k
            gfid  = int((LL-1)/(dim_sr*dim_sq*dim_a*dim_gq))+1
            LL    = modulo((LL-1),dim_sr*dim_sq*dim_a*dim_gq)+1
            gqid  = int((LL-1)/(dim_sr*dim_sq*dim_a))+1
            LL    = modulo((LL-1),dim_sr*dim_sq*dim_a)+1
            aid   = int((LL-1)/(dim_sr*dim_sq))+1
            LL    = modulo((LL-1),dim_sr*dim_sq)+1
            sqid  = int((LL-1)/(dim_sr))+1
            LL    = modulo((LL-1),dim_sr)+1
            srid  = LL

            vthc = sqrt(Ti0_1D(aid))
            ! mesh_cell_vol = (cell spatial volume)*integral(2*pi*sr^2*sin(sq))
            !               = (cell spatial volume)*2*pi*integral(sr^2)_(ur0-dsr/2)^(ur0+dsr/2) * integral(sin(sq))_(uq0-dsq/2)^(uq0+dsq/2)
            ! Cell_dVol_x3 = PtlSrcCV%Cell_dV(gqid,aid)
            Cell_dVol_x3 = Cell_dV(gqid,aid)

            sq0 = max(sq_min, sq(sqid)-0.5*dsq)
            sq1 = min(sq_max, sq(sqid)+0.5*dsq)
            Cell_dVol_sq = (-cos(sq1)+cos(sq0))

            sr0 = max(sr_min, sr(srid)-0.5*dsr)
            sr1 = min(sr_max, sr(srid)+0.5*dsr)
            ptl_ur_bound(1) = tgt_ur_bound(aid,gqid,gfid,1)
            ptl_ur_bound(2) = tgt_ur_bound(aid,gqid,gfid,2)
            ur0 = ptl_ur_bound(1) + sr0*(ptl_ur_bound(2)-ptl_ur_bound(1))
            ur1 = ptl_ur_bound(1) + sr1*(ptl_ur_bound(2)-ptl_ur_bound(1))
            Cell_dVol_sr = (ur1**3.0-ur0**3.0)/3.0

            mesh_cell_vol = 2.0*pi*(vthc**3.0)*Cell_dVol_x3*Cell_dVol_sr*Cell_dVol_sq

            ! mesh_cell_area is without Jacobian by v_perp (2*pi*vperp) uspace
            mesh_cell_area =  (min(a_max,a_1D(aid)+0.5*da)-max(a_min,a_1D(aid)-0.5*da))*dgq*dgf*0.5*(ur1**2.0-ur0**2.0)*(sq1-sq0)

            ! Normalize using Cell_dxdv
            mesh_norm_dxdv(srid,sqid,aid,gqid,gfid) = mesh_norm_dxdv(srid,sqid,aid,gqid,gfid)/mesh_cell_vol
            mesh_MarkerDen(srid,sqid,aid,gqid,gfid) = mesh_MarkerDen(srid,sqid,aid,gqid,gfid)/mesh_cell_area
            cell_vol_5D(srid,sqid,aid,gqid,gfid) = mesh_cell_vol
            cell_area_5D(srid,sqid,aid,gqid,gfid) = mesh_cell_area
         enddo
         deallocate(Ti0_1D,Cell_dV)

         !$acc update self(mesh_norm_dxdv,mesh_MarkerDen)

         allocate(tmp_mesh(dim_sr,dim_sq,dim_a,dim_gq,dim_gf))
         ! Communicate with other processors
         if(npartdom>1) then
            tmp_mesh = mesh_norm_dxdv
            length = dim_sr*dim_sq*dim_a*dim_gq*dim_gf
            call MPI_ALLREDUCE(tmp_mesh, mesh_norm_dxdv, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

            tmp_mesh = mesh_MarkerDen
            length = dim_sr*dim_sq*dim_a*dim_gq*dim_gf
            call MPI_ALLREDUCE(tmp_mesh, mesh_MarkerDen, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
         endif

         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically although they have different gqid)
         mesh_norm_dxdv(:,:,:,1,:) = mesh_norm_dxdv(:,:,:,1,:) + mesh_norm_dxdv(:,:,:,dim_gq,:)
         mesh_norm_dxdv(:,:,:,dim_gq,:) = mesh_norm_dxdv(:,:,:,1,:)

         mesh_MarkerDen(:,:,:,1,:) = mesh_MarkerDen(:,:,:,1,:) + mesh_MarkerDen(:,:,:,dim_gq,:)
         mesh_MarkerDen(:,:,:,dim_gq,:) = mesh_MarkerDen(:,:,:,1,:)

         !toroidal communication
         length = dim_sr*dim_sq*dim_a*dim_gq
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1, tmp_mesh(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         mesh_norm_dxdv(:,:,:,:,dim_gf) = mesh_norm_dxdv(:,:,:,:,dim_gf) + tmp_mesh(:,:,:,:,dim_gf)
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1, mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

         call MPI_SENDRECV(mesh_MarkerDen(:,:,:,:,1),length,mpi_Rsize,left_pe,1, tmp_mesh(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         mesh_MarkerDen(:,:,:,:,dim_gf) = mesh_MarkerDen(:,:,:,:,dim_gf) + tmp_mesh(:,:,:,:,dim_gf)
         call MPI_SENDRECV(mesh_MarkerDen(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1, mesh_MarkerDen(:,:,:,:,1),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)



         !Write into hdf5 file
         if(mype==0 .and. flag_debug_output==.true.) then
            write(filename, '("mesh_norm_Ion_Adap_",i5.5,"_a=",F0.1,".h5")') istep/msnap, 0.5*(a_min+a_max)

            call h5open_f(h5err)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err)
            if(h5err==-1) then ! If there is no file
               call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)
            else ! If there is file
            write(filename, '("mesh_norm_Ion_Adap_after_",i5.5,"_a=",F0.1,".h5")') istep/msnap, 0.5*(a_min+a_max)
               call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)
            endif

            gid=file_id

            h5dims(1) = dim_sr
            call h5ltmake_dataset_f(gid, "sr",1,h5dims,H5T_NATIVE_DOUBLE, sr, h5err)
            h5dims(1) = dim_sq
            call h5ltmake_dataset_f(gid, "sq",1,h5dims,H5T_NATIVE_DOUBLE, sq, h5err)
            h5dims(1) = dim_a
            call h5ltmake_dataset_f(gid, "a",1,h5dims,H5T_NATIVE_DOUBLE, a_1D, h5err)
            h5dims(1) = dim_gq
            call h5ltmake_dataset_f(gid, "gq",1,h5dims,H5T_NATIVE_DOUBLE, gq, h5err)
            h5dims(1) = dim_gf
            call h5ltmake_dataset_f(gid, "gf",1,h5dims,H5T_NATIVE_DOUBLE, gf, h5err)


            allocate(ur_4D(dim_sr,dim_a,dim_gq,dim_gf))
            allocate(spara_2D(dim_sr,dim_sq),sperp_2D(dim_sr,dim_sq))
            allocate(upara_5D(dim_sr,dim_sq,dim_a,dim_gq,dim_gf),uperp_5D(dim_sr,dim_sq,dim_a,dim_gq,dim_gf))

            do gfid=1,dim_gf
               do gqid=1,dim_gq
                  do aid=1,dim_a
                     do srid=1,dim_sr
                        ur_4D(srid,aid,gqid,gfid) = tgt_ur_bound(aid,gqid,gfid,1)+ sr(srid)*(tgt_ur_bound(aid,gqid,gfid,2)-tgt_ur_bound(aid,gqid,gfid,1))
                        do sqid=1,dim_sq
                           upara_5D(srid,sqid,aid,gqid,gfid) = ur_4D(srid,aid,gqid,gfid)*cos(sq(sqid))
                           uperp_5D(srid,sqid,aid,gqid,gfid) = ur_4D(srid,aid,gqid,gfid)*sin(sq(sqid))
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            do srid=1,dim_sr
               do sqid=1,dim_sq
                  spara_2D(srid,sqid) = sr(srid)*cos(sq(sqid))
                  sperp_2D(srid,sqid) = sr(srid)*sin(sq(sqid))
               enddo
            enddo

            h5dims(1:2) = [dim_sr,dim_sq]
            call h5ltmake_dataset_f(gid, "spara_2D",2,h5dims,H5T_NATIVE_DOUBLE, spara_2D, h5err)
            call h5ltmake_dataset_f(gid, "sperp_2D",2,h5dims,H5T_NATIVE_DOUBLE, sperp_2D, h5err)

            h5dims(1:4) = [dim_sr,dim_a,dim_gq,dim_gf]
            call h5ltmake_dataset_f(gid, "ur_4D",4,h5dims,H5T_NATIVE_DOUBLE, ur_4D, h5err)

            h5dims(1:5) = [dim_sr,dim_sq,dim_a,dim_gq,dim_gf]
            call h5ltmake_dataset_f(gid, "mesh_norm_dxdv",5,h5dims,H5T_NATIVE_DOUBLE, mesh_norm_dxdv, h5err)
            call h5ltmake_dataset_f(gid, "mesh_MarkerDen",5,h5dims,H5T_NATIVE_DOUBLE, mesh_MarkerDen, h5err)
            call h5ltmake_dataset_f(gid, "mesh_cell_vol",5,h5dims,H5T_NATIVE_DOUBLE, cell_vol_5D, h5err)
            call h5ltmake_dataset_f(gid, "mesh_cell_area",5,h5dims,H5T_NATIVE_DOUBLE, cell_area_5D, h5err)
            call h5ltmake_dataset_f(gid, "upara_5D",5,h5dims,H5T_NATIVE_DOUBLE, upara_5D, h5err)
            call h5ltmake_dataset_f(gid, "uperp_5D",5,h5dims,H5T_NATIVE_DOUBLE, uperp_5D, h5err)


            tmp_mesh = mesh_norm_dxdv
            call Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions([10,10],tmp_mesh)
            call h5ltmake_dataset_f(gid, "mesh_norm_dxdv_SM",5,h5dims,H5T_NATIVE_DOUBLE, tmp_mesh, h5err)
            tmp_mesh = mesh_MarkerDen
            call Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions([10,10],tmp_mesh)
            call h5ltmake_dataset_f(gid, "mesh_MarkerDen_SM",5,h5dims,H5T_NATIVE_DOUBLE, tmp_mesh, h5err)


            call h5gclose_f(gid, h5err)
            call h5fclose_f(file_id, h5err)
            call h5close_f(h5err)

            deallocate(spara_2D,sperp_2D)
            deallocate(ur_4D)
            deallocate(upara_5D,uperp_5D)
         endif

         deallocate(a_1D,gq,gf)
         deallocate(sr,sq)
         deallocate(tmp_mesh)
         deallocate(cell_vol_5D)
         deallocate(cell_area_5D)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Global_mesh_norm_dxdv_of_Ions_PolarVelCoord_on_TgtUr takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Measure_Global_mesh_norm_dxdv_of_Ions_PolarVelCoord_on_TgtUr


      subroutine Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord_on_TgtUr(dim_sr, dim_sq, tgt_ur_bound, mesh_norm_dxdv, mesh_MarkerDen, InFlag_debug_output)
         use global_parameters
         use particle_array
         use particle_decomp
         use Structured_Grid
         use esi_interface
         use mpi_defs
         use function_gg_module
         use hdf5
         use h5lt
         implicit none
         integer, intent(in) :: dim_sr, dim_sq
         real(wp), intent(in) :: tgt_ur_bound(:,:,:,:) ! (dim_a,dim_gq,dim_gf,2); upper ur and lower ur bound on 3-D grid (a,gq,gf)
         real(wp), intent(inout) :: mesh_norm_dxdv(:,:,:,:,:) ! 5-D grid (sr,sq,a,gq,gf) normalized dxdv by pre-occupying particles
         real(wp), intent(inout) :: mesh_MarkerDen(:,:,:,:,:) ! 5-D grid (sr,sq,a,gq,gf) normalized density by pre-occupying particles (on upace, no Jacobian)
         logical, intent(in), optional :: InFlag_debug_output
         logical :: flag_debug_output
         !$acc declare create(mesh_norm_dxdv,mesh_MarkerDen)
         real(wp), allocatable :: tmp_mesh(:,:,:,:,:),cell_vol_5D(:,:,:,:,:), cell_area_5D(:,:,:,:,:)
         real(wp), allocatable :: sr(:), sq(:)
         real(wp), allocatable :: a_1D(:), gq(:), gf(:)
         integer :: tmp_aid, tgt_aid, tgt_gqid, tgt_gfid, tgt_srid, tgt_sqid
         real(wp) ptl_ur, ptl_sr, ptl_sq, ptl_ur_bound(2)
         real(wp) tgt_gq, ws0,ws1, P_safety_q
         real(wp) upara,uperp, Cell_dVol_x3, Cell_dVol_sr, Cell_dVol_sq
         real(wp), allocatable :: spara_2D(:,:), sperp_2D(:,:), ur_4D(:,:,:,:), upara_5D(:,:,:,:,:), uperp_5D(:,:,:,:,:)
         integer m,i,j,k,LL,length,ierror,istatus(MPI_STATUS_SIZE)
         integer aid,gqid,gfid,srid,sqid
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max, sr_min, sr_max, sq_min, sq_max
         real(wp) :: da, dgq, dgf, dsr, dsq, inv_da, inv_dgq, inv_dgf, inv_dsr, inv_dsq
         real(wp) :: sr0,sr1,sq0,sq1,ur0,ur1
         integer :: dim_a, dim_gq, dim_gf
         real(wp), external :: temperature_e
         real(wp) mesh_cell_vol, mesh_cell_area ! 5-D mesh volume and area of each mesh grid
         real(wp) time0,time1
         real(wp), allocatable :: Te0_1D(:), Cell_dV(:,:)

         character(len=128) filename, gname
         integer(HID_T) file_id, gid
         integer :: h5err
         integer(HSIZE_T) :: h5dims(5)

         ! Global mesh means that the mesh covers all the spatial range (a0-a1, 0-2pi, 0-2pi)
         ! by taking into account all processors and their communications.
         ! Also, the periodical conditions into theta and toroidal direcitons are considered

#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         flag_debug_output = .false.
         if(present(InFlag_debug_output)) flag_debug_output=InFlag_debug_output

         ! ! For simplicity, let's use a zero-th order shape function (closest near grid)
         if(acc_newer_H_zeon) then
            !$acc update device(zeon(:,:me),zeon0(:,:me),RK_cum_zeon(:,:me))
            acc_newer_H_zeon=.false.
            acc_newer_D_zeon=.false.
         endif

         call n_T_u_group_e   ! get Te, note that back_f0_e cannot be calculated here because it requires zeon(4) and zeon(6)
         call esigetb(zeon,me,zeondim)   ! get B-->BB



         ! allocate a,gq,gf and copy SG information into
         allocate(a_1D(PSCV%dim_a),gq(PSCV%dim_gq),gf(2))
         allocate(sr(dim_sr),sq(dim_sq))
         a_1D=PSCV%a_1D
         gq=PSCV%gq_1D
         gf(1)=zetamin
         gf(2)=zetamax

         sr_min = 0; sr_max = 1;
         sq_min = 0; sq_max = pi;
         dsr=(sr_max-sr_min)/(dim_sr-1); inv_dsr = 1.0/dsr
         dsq=(sq_max-sq_min)/(dim_sq-1); inv_dsq = 1.0/dsq
         do i=1,dim_sr
            sr(i) = sr_min + (i-1)*dsr
         enddo
         do i=1,dim_sq
            sq(i) = sq_min + (i-1)*dsq
         enddo


         a_min = minval(a_1D); a_max = maxval(a_1D); dim_a = size(a_1D)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         da=(a_max-a_min)/(dim_a-1);     inv_da  = 1.0/da
         dgq=(gq_max-gq_min)/(dim_gq-1); inv_dgq = 1.0/dgq
         dgf=(gf_max-gf_min)/(dim_gf-1); inv_dgf = 1.0/dgf

         ! Measure normalized phase volume from electrons that reside in mesh domain
         ! mesh_norm_dxdv was allocated at the outside of this subroutine
         mesh_norm_dxdv=0.0
         mesh_MarkerDen=0.0
         !$acc update device(mesh_norm_dxdv,mesh_MarkerDen)

         !$acc parallel loop present(mesh_norm_dxdv,mesh_MarkerDen) private(ptl_ur_bound)
         do m=1,me

            if( (zeon(1,m)-a_min)*(zeon(1,m)-a_max)>0 ) then
               ! If the particle is outside the PSCV domain, cycle the loop
               cycle
            endif

            tgt_aid  = ceiling((zeon(1,m)-a_min)*inv_da+0.5) ! Find nearst grid
            tgt_gfid = ceiling((zeon(3,m)-gf_min)*inv_dgf+0.5) ! Find nearst grid

            ! Interpolate safety factor
            tmp_aid = max(1,min(int((zeon(1,m)-PSCV%a_min)*PSCV%inv_da) +1,PSCV%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
            ws1 = (zeon(1,m)-PSCV%a_1D(tmp_aid))*PSCV%inv_da
            ws0 = (1.0-ws1)
            P_safety_q = ws0*PSCV%safety_q_1D(tmp_aid) + ws1*PSCV%safety_q_1D(tmp_aid+1)

            tgt_gq = modulo(zeon(2,m) + (gf(tgt_gfid)-zeon(3,m))/P_safety_q, 2.0*pi)
            tgt_gqid = ceiling((tgt_gq-gq_min)*inv_dgq+0.5) ! Find nearst grid

            ! Find ur bound of the particle
            ptl_ur_bound(1) = tgt_ur_bound(tgt_aid,tgt_gqid,tgt_gfid,1)
            ptl_ur_bound(2) = tgt_ur_bound(tgt_aid,tgt_gqid,tgt_gfid,2)

            ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
            ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))
            upara = zeon(4,m)*csc0*BB(m)*sqrt(m0ms/(msm0*backTe(m)))! upara
            uperp = zeon(6,m)*sqrt(2.0*BB(m))/sqrt(msm0*backTe(m))  ! uperp

            ptl_ur = sqrt(upara**2.0+uperp**2.0)

            if( (ptl_ur-ptl_ur_bound(1))*(ptl_ur-ptl_ur_bound(2)) <0) then
               ! Measure particle's contribution to the norm_phase_volume and MarkerDen
               ptl_sr = (ptl_ur-ptl_ur_bound(1))/(ptl_ur_bound(2)-ptl_ur_bound(1))
               tgt_srid = ceiling((ptl_sr-sr_min)*inv_dsr+0.5)
               ptl_sq = acos(upara/ptl_ur)
               tgt_sqid = ceiling((ptl_sq-sq_min)*inv_dsq+0.5)

               !$acc atomic update
               mesh_norm_dxdv(tgt_srid,tgt_sqid,tgt_aid,tgt_gqid,tgt_gfid) = mesh_norm_dxdv(tgt_srid,tgt_sqid,tgt_aid,tgt_gqid,tgt_gfid) + zeon0(6,m)
               !$acc atomic update
               mesh_MarkerDen(tgt_srid,tgt_sqid,tgt_aid,tgt_gqid,tgt_gfid) = mesh_MarkerDen(tgt_srid,tgt_sqid,tgt_aid,tgt_gqid,tgt_gfid) + 1.0
            endif
         enddo


         allocate(Te0_1D(PSCV%dim_a),Cell_dV(PSCV%dim_gq,PSCV%dim_a))

         Te0_1D = PSCV%Te0_1D
         Cell_dV = PSCV%Cell_dV



         allocate(cell_vol_5D(dim_sr,dim_sq,dim_a,dim_gq,dim_gf))
         allocate(cell_area_5D(dim_sr,dim_sq,dim_a,dim_gq,dim_gf))

         length=dim_sq*dim_sr*dim_a*dim_gq*dim_gf
         ! Normalize the results
         !$acc parallel loop gang vector present(mesh_norm_dxdv,mesh_MarkerDen) private(LL,k,ptl_ur_bound)
         do k=1,length
            LL=k
            gfid  = int((LL-1)/(dim_sr*dim_sq*dim_a*dim_gq))+1
            LL    = modulo((LL-1),dim_sr*dim_sq*dim_a*dim_gq)+1
            gqid  = int((LL-1)/(dim_sr*dim_sq*dim_a))+1
            LL    = modulo((LL-1),dim_sr*dim_sq*dim_a)+1
            aid   = int((LL-1)/(dim_sr*dim_sq))+1
            LL    = modulo((LL-1),dim_sr*dim_sq)+1
            sqid  = int((LL-1)/(dim_sr))+1
            LL    = modulo((LL-1),dim_sr)+1
            srid  = LL

            vthc = sqrt(Te0_1D(aid)*m0ms)
               ! mesh_cell_vol = (cell spatial volume)*integral(2*pi*sr^2*sin(sq))
               !               = (cell spatial volume)*2*pi*integral(sr^2)_(ur0-dsr/2)^(ur0+dsr/2) * integral(sin(sq))_(uq0-dsq/2)^(uq0+dsq/2)
            ! Cell_dVol_x3 = PtlSrcCV%Cell_dV(gqid,aid)
            Cell_dVol_x3 = Cell_dV(gqid,aid)


            sq0 = max(sq_min, sq(sqid)-0.5*dsq)
            sq1 = min(sq_max, sq(sqid)+0.5*dsq)
            Cell_dVol_sq = (-cos(sq1)+cos(sq0))


            sr0 = max(sr_min, sr(srid)-0.5*dsr)
            sr1 = min(sr_max, sr(srid)+0.5*dsr)
            ptl_ur_bound(1) = tgt_ur_bound(aid,gqid,gfid,1)
            ptl_ur_bound(2) = tgt_ur_bound(aid,gqid,gfid,2)
            ur0 = ptl_ur_bound(1) + sr0*(ptl_ur_bound(2)-ptl_ur_bound(1))
            ur1 = ptl_ur_bound(1) + sr1*(ptl_ur_bound(2)-ptl_ur_bound(1))
            Cell_dVol_sr = (ur1**3.0-ur0**3.0)/3.0

            mesh_cell_vol = 2.0*pi*(vthc**3.0)*Cell_dVol_x3*Cell_dVol_sr*Cell_dVol_sq

            ! mesh_cell_area is without Jacobian by v_perp (2*pi*vperp) uspace
            mesh_cell_area =  (min(a_max,a_1D(aid)+0.5*da)-max(a_min,a_1D(aid)-0.5*da))*dgq*dgf*0.5*(ur1**2.0-ur0**2.0)*(sq1-sq0)

            ! Normalize using Cell_dxdv
            mesh_norm_dxdv(srid,sqid,aid,gqid,gfid) = mesh_norm_dxdv(srid,sqid,aid,gqid,gfid)/mesh_cell_vol
            mesh_MarkerDen(srid,sqid,aid,gqid,gfid) = mesh_MarkerDen(srid,sqid,aid,gqid,gfid)/mesh_cell_area
            cell_vol_5D(srid,sqid,aid,gqid,gfid) = mesh_cell_vol
            cell_area_5D(srid,sqid,aid,gqid,gfid) = mesh_cell_area
         enddo
         deallocate(Te0_1D,Cell_dV)

         !$acc update self(mesh_norm_dxdv,mesh_MarkerDen)

         allocate(tmp_mesh(dim_sr,dim_sq,dim_a,dim_gq,dim_gf))
         ! Communicate with other processors
         if(npartdom>1) then
            tmp_mesh = mesh_norm_dxdv
            length = dim_sr*dim_sq*dim_a*dim_gq*dim_gf
            call MPI_ALLREDUCE(tmp_mesh, mesh_norm_dxdv, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)

            tmp_mesh = mesh_MarkerDen
            length = dim_sr*dim_sq*dim_a*dim_gq*dim_gf
            call MPI_ALLREDUCE(tmp_mesh, mesh_MarkerDen, length, mpi_Rsize, MPI_SUM, partd_comm, ierror)
         endif

         ! Match poloidal end cell
         ! (gq=0 and gq=2pi is the same position physically although they have different gqid)
         mesh_norm_dxdv(:,:,:,1,:) = mesh_norm_dxdv(:,:,:,1,:) + mesh_norm_dxdv(:,:,:,dim_gq,:)
         mesh_norm_dxdv(:,:,:,dim_gq,:) = mesh_norm_dxdv(:,:,:,1,:)

         mesh_MarkerDen(:,:,:,1,:) = mesh_MarkerDen(:,:,:,1,:) + mesh_MarkerDen(:,:,:,dim_gq,:)
         mesh_MarkerDen(:,:,:,dim_gq,:) = mesh_MarkerDen(:,:,:,1,:)

         !toroidal communication
         length = dim_sr*dim_sq*dim_a*dim_gq
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1, tmp_mesh(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         mesh_norm_dxdv(:,:,:,:,dim_gf) = mesh_norm_dxdv(:,:,:,:,dim_gf) + tmp_mesh(:,:,:,:,dim_gf)
         call MPI_SENDRECV(mesh_norm_dxdv(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1, mesh_norm_dxdv(:,:,:,:,1),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)

         call MPI_SENDRECV(mesh_MarkerDen(:,:,:,:,1),length,mpi_Rsize,left_pe,1, tmp_mesh(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1,toroidal_comm,istatus,ierror)
         mesh_MarkerDen(:,:,:,:,dim_gf) = mesh_MarkerDen(:,:,:,:,dim_gf) + tmp_mesh(:,:,:,:,dim_gf)
         call MPI_SENDRECV(mesh_MarkerDen(:,:,:,:,dim_gf),length,mpi_Rsize,right_pe,1, mesh_MarkerDen(:,:,:,:,1),length,mpi_Rsize,left_pe,1,toroidal_comm,istatus,ierror)



         !Write into hdf5 file
         if(mype==0 .and. flag_debug_output==.true.) then
            ! write(filename, '("polar_mesh_norm_Elec_dxdv_",i5.5,"_urmin=",F0.2,".h5")') istep/msnap, sr_min
            ! write(filename, '("mesh_norm_Elec_Adap_",i5.5,".h5")') istep/msnap
            write(filename, '("mesh_norm_Elec_Adap_",i5.5,"_a=",F0.1,".h5")') istep/msnap, 0.5*(a_min+a_max)

            call h5open_f(h5err)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err)
            if(h5err==-1) then ! If there is no file
               call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)
            else ! If there is file
               ! write(filename, '("mesh_norm_Elec_Adap_after",i5.5,".h5")') istep/msnap
            write(filename, '("mesh_norm_Elec_Adap_after_",i5.5,"_a=",F0.1,".h5")') istep/msnap, 0.5*(a_min+a_max)
               call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)
            endif

            ! write(gname,'("polar")')
            ! call h5gcreate_f(file_id, gname, gid, h5err)
            gid=file_id


            h5dims(1) = dim_sr
            call h5ltmake_dataset_f(gid, "sr",1,h5dims,H5T_NATIVE_DOUBLE, sr, h5err)
            h5dims(1) = dim_sq
            call h5ltmake_dataset_f(gid, "sq",1,h5dims,H5T_NATIVE_DOUBLE, sq, h5err)
            h5dims(1) = dim_a
            call h5ltmake_dataset_f(gid, "a",1,h5dims,H5T_NATIVE_DOUBLE, a_1D, h5err)
            h5dims(1) = dim_gq
            call h5ltmake_dataset_f(gid, "gq",1,h5dims,H5T_NATIVE_DOUBLE, gq, h5err)
            h5dims(1) = dim_gf
            call h5ltmake_dataset_f(gid, "gf",1,h5dims,H5T_NATIVE_DOUBLE, gf, h5err)


            allocate(ur_4D(dim_sr,dim_a,dim_gq,dim_gf))
            allocate(spara_2D(dim_sr,dim_sq),sperp_2D(dim_sr,dim_sq))
            allocate(upara_5D(dim_sr,dim_sq,dim_a,dim_gq,dim_gf),uperp_5D(dim_sr,dim_sq,dim_a,dim_gq,dim_gf))

            do gfid=1,dim_gf
               do gqid=1,dim_gq
                  do aid=1,dim_a
                     do srid=1,dim_sr
                        ur_4D(srid,aid,gqid,gfid) = tgt_ur_bound(aid,gqid,gfid,1)+ sr(srid)*(tgt_ur_bound(aid,gqid,gfid,2)-tgt_ur_bound(aid,gqid,gfid,1))
                        do sqid=1,dim_sq
                           upara_5D(srid,sqid,aid,gqid,gfid) = ur_4D(srid,aid,gqid,gfid)*cos(sq(sqid))
                           uperp_5D(srid,sqid,aid,gqid,gfid) = ur_4D(srid,aid,gqid,gfid)*sin(sq(sqid))
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            do srid=1,dim_sr
               do sqid=1,dim_sq
                  spara_2D(srid,sqid) = sr(srid)*cos(sq(sqid))
                  sperp_2D(srid,sqid) = sr(srid)*sin(sq(sqid))
               enddo
            enddo

            h5dims(1:2) = [dim_sr,dim_sq]
            call h5ltmake_dataset_f(gid, "spara_2D",2,h5dims,H5T_NATIVE_DOUBLE, spara_2D, h5err)
            call h5ltmake_dataset_f(gid, "sperp_2D",2,h5dims,H5T_NATIVE_DOUBLE, sperp_2D, h5err)

            h5dims(1:4) = [dim_sr,dim_a,dim_gq,dim_gf]
            call h5ltmake_dataset_f(gid, "ur_4D",4,h5dims,H5T_NATIVE_DOUBLE, ur_4D, h5err)

            h5dims(1:5) = [dim_sr,dim_sq,dim_a,dim_gq,dim_gf]
            call h5ltmake_dataset_f(gid, "mesh_norm_dxdv",5,h5dims,H5T_NATIVE_DOUBLE, mesh_norm_dxdv, h5err)
            call h5ltmake_dataset_f(gid, "mesh_MarkerDen",5,h5dims,H5T_NATIVE_DOUBLE, mesh_MarkerDen, h5err)
            call h5ltmake_dataset_f(gid, "mesh_cell_vol",5,h5dims,H5T_NATIVE_DOUBLE, cell_vol_5D, h5err)
            call h5ltmake_dataset_f(gid, "mesh_cell_area",5,h5dims,H5T_NATIVE_DOUBLE, cell_area_5D, h5err)
            call h5ltmake_dataset_f(gid, "upara_5D",5,h5dims,H5T_NATIVE_DOUBLE, upara_5D, h5err)
            call h5ltmake_dataset_f(gid, "uperp_5D",5,h5dims,H5T_NATIVE_DOUBLE, uperp_5D, h5err)


            tmp_mesh = mesh_norm_dxdv
            call Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions([10,10],tmp_mesh)
            call h5ltmake_dataset_f(gid, "mesh_norm_dxdv_SM",5,h5dims,H5T_NATIVE_DOUBLE, tmp_mesh, h5err)
            tmp_mesh = mesh_MarkerDen
            call Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions([10,10],tmp_mesh)
            call h5ltmake_dataset_f(gid, "mesh_MarkerDen_SM",5,h5dims,H5T_NATIVE_DOUBLE, tmp_mesh, h5err)


            call h5gclose_f(gid, h5err)
            call h5fclose_f(file_id, h5err)
            call h5close_f(h5err)

            deallocate(spara_2D,sperp_2D)
            deallocate(ur_4D)
            deallocate(upara_5D,uperp_5D)
         endif

         deallocate(a_1D,gq,gf)
         deallocate(sr,sq)
         deallocate(tmp_mesh)
         deallocate(cell_vol_5D)
         deallocate(cell_area_5D)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord_on_TgtUr takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Measure_Global_mesh_norm_dxdv_of_Electrons_PolarVelCoord_on_TgtUr


      subroutine Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions(num_smooth_gq_a,tgt_mesh5D, Reference_Vol_2D_gq_a)
         use global_parameters
         implicit none

         integer, intent(in) :: num_smooth_gq_a(:)
         real(wp),intent(inout) :: tgt_mesh5D(:,:,:,:,:)
         real(wp),optional:: Reference_Vol_2D_gq_a(PSCV%dim_gq, PSCV%dim_a)
         real(wp) :: VV(PSCV%dim_gq, PSCV%dim_a)
         integer :: smid,aid,qid,Na,Nq,Nsr,Nsq
         real(wp), allocatable :: SM_Data5D(:,:,:,:,:),ori_SM_Data5D(:,:,:,:,:) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
         real(wp) :: SM_gq_Coeff(-1:1,0:PSCV%dim_a+1, 0:PSCV%dim_gq+1), SM_a_Coeff(-1:1,0:PSCV%dim_a+1, 0:PSCV%dim_gq+1)

         Na = PSCV%dim_a
         Nq = PSCV%dim_gq

         Nsr = size(tgt_mesh5D,1)
         Nsq = size(tgt_mesh5D,2)

         if(Na/=size(tgt_mesh5D,3) .or. Nq/=size(tgt_mesh5D,4)) then
            write(stdout,'((A,5I0,A))') "Error! tgt_mesh5D does not fit the PtlSrcCV!! dims(tgt_mesh5D)=(", size(tgt_mesh5D,1), size(tgt_mesh5D,2), size(tgt_mesh5D,3), size(tgt_mesh5D,4), size(tgt_mesh5D,5),")"
            stop
         endif

         allocate(SM_Data5D(Nsr,Nsq,0:Na+1,0:Nq+1,2), ori_SM_Data5D(Nsr,Nsq,0:Na+1,0:Nq+1,2))

         VV=1.0
         if(present(Reference_Vol_2D_gq_a)) VV=Reference_Vol_2D_gq_a

         ! Get Smooth Coeffs which conserves Volume Weighted Integration
         call Get_Conserving_Smooth_Coeffs_on_PtlSrcCV(VV, SM_a_Coeff, SM_gq_Coeff)


         if(size(num_smooth_gq_a)==1) then
            SM_Data5D = 0.0
            SM_Data5D(:,:,1:Na,1:Nq,:) = tgt_mesh5D
            do smid=1,num_smooth_gq_a(1)
               ! Set ghost cell periodic in theta direction
               SM_Data5D(:,:,:,0,:)    = SM_Data5D(:,:,:,Nq-1,:)
               SM_Data5D(:,:,:,Nq+1,:) = SM_Data5D(:,:,:,2,:)
               ! Smooth in theta direction
               ori_SM_Data5D = SM_Data5D
               do aid=1,Na
                  do qid=1,Nq
                     SM_Data5D(:,:,aid,qid,:) = SM_gq_Coeff(-1,aid,qid)*ori_SM_Data5D(:,:,aid,qid-1,:) + SM_gq_Coeff(0,aid,qid)*ori_SM_Data5D(:,:,aid,qid,:) + SM_gq_Coeff(1,aid,qid)*ori_SM_Data5D(:,:,aid,qid+1,:)
                  enddo
               enddo
               ! Smooth in radial direction
               ori_SM_Data5D = SM_Data5D
               do aid=1,Na
                  do qid=1,Nq
                     SM_Data5D(:,:,aid,qid,:) = SM_a_Coeff(-1,aid,qid)*ori_SM_Data5D(:,:,aid-1,qid,:) + SM_a_Coeff(0,aid,qid)*ori_SM_Data5D(:,:,aid,qid,:) + SM_a_Coeff(1,aid,qid)*ori_SM_Data5D(:,:,aid+1,qid,:)
                  enddo
               enddo
            enddo
            tgt_mesh5D = SM_Data5D(:,:,1:Na,1:Nq,:)
         elseif(size(num_smooth_gq_a)==2) then
            SM_Data5D = 0.0
            SM_Data5D(:,:,1:Na,1:Nq,:) = tgt_mesh5D
            do smid=1,num_smooth_gq_a(1)
               ! Set ghost cell periodic in theta direction
               SM_Data5D(:,:,:,0,:)    = SM_Data5D(:,:,:,Nq-1,:)
               SM_Data5D(:,:,:,Nq+1,:) = SM_Data5D(:,:,:,2,:)
               ! Smooth in theta direction
               ori_SM_Data5D = SM_Data5D
               do aid=1,Na
                  do qid=1,Nq
                     SM_Data5D(:,:,aid,qid,:) = SM_gq_Coeff(-1,aid,qid)*ori_SM_Data5D(:,:,aid,qid-1,:) + SM_gq_Coeff(0,aid,qid)*ori_SM_Data5D(:,:,aid,qid,:) + SM_gq_Coeff(1,aid,qid)*ori_SM_Data5D(:,:,aid,qid+1,:)
                     ! SM_Data5D(:,:,aid,qid,:) = 0.25*ori_SM_Data5D(:,:,aid,qid-1,:) + 0.5*ori_SM_Data5D(:,:,aid,qid,:) + 0.25*ori_SM_Data5D(:,:,aid,qid+1,:)
                  enddo
               enddo
            enddo
            do smid=1,num_smooth_gq_a(2)
               ! Smooth in radial direction
               ori_SM_Data5D = SM_Data5D
               do aid=1,Na
                  do qid=1,Nq
                     SM_Data5D(:,:,aid,qid,:) = SM_a_Coeff(-1,aid,qid)*ori_SM_Data5D(:,:,aid-1,qid,:) + SM_a_Coeff(0,aid,qid)*ori_SM_Data5D(:,:,aid,qid,:) + SM_a_Coeff(1,aid,qid)*ori_SM_Data5D(:,:,aid+1,qid,:)
                     ! SM_Data5D(:,:,aid,qid,:) = 0.25*ori_SM_Data5D(:,:,aid-1,qid,:) + 0.5*ori_SM_Data5D(:,:,aid,qid,:) + 0.25*ori_SM_Data5D(:,:,aid+1,qid,:)
                  enddo
               enddo
            enddo
            tgt_mesh5D = SM_Data5D(:,:,1:Na,1:Nq,:)
         else
            write(stdout,*) "Error! Check num_smooth_gq_a = ",num_smooth_gq_a(:)
            stop
         endif

         deallocate(SM_Data5D, ori_SM_Data5D)

      end subroutine Smooth_PtlSrcCV_mesh5D_along_radial_and_poloidal_directions

      subroutine Get_Conserving_Smooth_Coeffs_on_PtlSrcCV(in_VOLorAREA, SM_a_Coeff, SM_gq_Coeff)
         implicit none

         integer :: Na,Nq
         real(wp), intent(in) :: in_VOLorAREA(PSCV%dim_a,PSCV%dim_gq) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
         real(wp) :: VV(0:PSCV%dim_a+1, 0:PSCV%dim_gq+1) ! Note that there are ghost cells surrouding boundaries (index 0, dim_*+1)
         real(wp), intent(out):: SM_gq_Coeff(-1:1, 0:PSCV%dim_a+1,0:PSCV%dim_gq+1), SM_a_Coeff(-1:1, 0:PSCV%dim_a+1,0:PSCV%dim_gq+1)

         Na = PSCV%dim_a
         Nq = PSCV%dim_gq

         VV = 0
         VV(1:Na,1:Nq) = in_VOLorAREA

         ! Note that theta direction is periodic
         VV(:,0) = VV(:,Nq-1)
         VV(:,Nq+1) = VV(:,2)

         SM_gq_Coeff(-1,:,1:Nq) = 0.5*VV(:,0:Nq-1)/(VV(:,0:Nq-1) + VV(:,1:Nq)) ! Coeff_L = 0.5*VL/(VL+VC)
         SM_gq_Coeff( 0,:,1:Nq) = 0.5*VV(:,1:Nq)*(VV(:,0:Nq-1)+2.0*VV(:,1:Nq)+VV(:,2:Nq+1))/((VV(:,0:Nq-1)+VV(:,1:Nq))*(VV(:,1:Nq)+VV(:,2:Nq+1))) ! 0.5*VC*(VL+2*VC+VR)/((VL+VC)*(VC+VR))
         SM_gq_Coeff( 1,:,1:Nq) = 0.5*VV(:,2:Nq+1)/(VV(:,2:Nq+1) + VV(:,1:Nq)) ! Coeff_R = 0.5*VR/(VR+VC)

         SM_a_Coeff(-1,1:Na,:) = 0.5*VV(0:Na-1,:)/(VV(0:Na-1,:) + VV(1:Na,:)) ! Coeff_L = 0.5*VL/(VL+VC)
         SM_a_Coeff( 0,1:Na,:) = 0.5*VV(1:Na,:)*(VV(0:Na-1,:)+2.0*VV(1:Na,:)+VV(2:Na+1,:))/((VV(0:Na-1,:)+VV(1:Na,:))*(VV(1:Na,:)+VV(2:Na+1,:))) ! 0.5*VC*(VL+2*VC+VR)/((VL+VC)*(VC+VR))
         SM_a_Coeff( 1,1:Na,:) = 0.5*VV(2:Na+1,:)/(VV(2:Na+1,:) + VV(1:Na,:)) ! Coeff_R = 0.5*VR/(VR+VC)

         ! To remove NsN values at the raidal ghost cells
         SM_gq_Coeff(:,0,:)   = 0
         SM_gq_Coeff(:,Na+1,:)= 0
         SM_a_Coeff(:,0,:)    = 0
         SM_a_Coeff(:,Na+1,:) = 0

      end subroutine Get_Conserving_Smooth_Coeffs_on_PtlSrcCV




      subroutine Generate_Ions_in_given_Phase_Space(upara,uperp,a,gq,gf,this_N_newP, mesh_norm_dxdv, in_label_val)
         use global_parameters
         use particle_array
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         use particle_decomp
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         real(wp), intent(in), dimension(:) :: upara,uperp,a,gq,gf
         integer, intent(in) :: this_N_newP
         real(wp), intent(in), optional :: mesh_norm_dxdv(:,:,:,:,:) ! 5-D grid (upara,uperp,a,gq,gf) normalized dxdv by pre-occupying particles
         real(wp), intent(in), optional :: in_label_val

         real(wp) label_val
         integer m,m0,i, tot_N_newP
         real(wp) :: vpara, vperp, integral_f0, target_w_tot
         real(wp), allocatable :: gi_val(:)
         real(wp), allocatable :: rand_5d(:,:)
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max, upara_min, upara_max, uperp_min, uperp_max
         real(wp) :: Del_a, Del_gq, Del_gf, Del_upara, Del_uperp
         real(wp) :: da, dgq, dgf, dupara, duperp
         integer :: dim_a, dim_gq, dim_gf, dim_upara, dim_uperp
         integer :: minpts, ifail=1, ierror
         real(wp) :: ac, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian
         real(wp) time0,time1

         real(wp) :: mesh_min_x(5),mesh_inv_dx(5)
         real(wp) :: ptl_data(5)
         integer :: mesh_idx(5), mesh_dims(5)
         real(wp), allocatable :: target_norm_phase_vol(:)
         !$acc declare create(target_norm_phase_vol)

         call MPI_ALLREDUCE(this_N_newP, tot_N_newP,1,MPI_INTEGER, MPI_SUM, partd_comm,ierror)

         if(tot_N_newP<=0) then
            return
         endif

         label_val=0
         if(present(in_label_val)) label_val = in_label_val

#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         upara_min = minval(upara); upara_max = maxval(upara); dim_upara = size(upara)
         uperp_min = minval(uperp); uperp_max = maxval(uperp); dim_uperp = size(uperp)
         a_min = minval(a); a_max = maxval(a); dim_a = size(a)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         Del_upara = upara_max-upara_min;  dupara=Del_upara/(dim_upara-1)
         Del_uperp = uperp_max-uperp_min;  duperp=Del_uperp/(dim_uperp-1)
         Del_a = a_max-a_min;       da=Del_a/(dim_a-1)
         Del_gq = gq_max-gq_min;    dgq=Del_gq/(dim_gq-1)
         Del_gf = gf_max-gf_min;    dgf=Del_gf/(dim_gf-1)


         ! Generate
         m0=mi+1
         mi=mi+this_N_newP

         allocate(gi_val(m0:mi), rand_5d(5,m0:mi))

         ! Prepare random numbers
         do m=m0,mi
            do i=1,5
               rand_5d(i,m) = sprng()
            enddo
         enddo


         do m=m0,mi
            zion(1,m) = a_min + rand_5d(1,m)*Del_a
            zion(7,m)=D01AHF(0.0_wp,zion(1,m),error_8,minpts,ac,dsda,-1,ifail)

            zion(2,m) = modulo(gq_min + rand_5d(2,m)*Del_gq,2.0*pi)  ! gq
            zion(3,m) = modulo(gf_min + rand_5d(3,m)*Del_gf,2.0*pi)  ! gf

            ilabels(m) = label_val
         enddo
         !$acc update device(zion(:,m0:mi),ilabels(m0:mi))

         call n_T_u_group(m0,mi)   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(backu(m0:mi),backT(m0:mi))
         call esigetbf_partial(zion,m0,mi,ziondim)
         backu(m0:mi)=backu(m0:mi)*F(m0:mi)/BB(m0:mi) ! =u_//=I*omega_t/B


         ! Uniform Loading
         do m=m0,mi
            vpara = (upara_min + rand_5d(4,m)*Del_upara)*sqrt(backT(m))
            vperp = (uperp_min + rand_5d(5,m)*Del_uperp)*sqrt(backT(m))
            gi_val(m) = tot_N_newP/(Del_a*Del_gq*Del_gf*Del_upara*Del_uperp*backT(m))/(Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)

            zion(4,m) = vpara/BB(m)
            zion(4,m) = zion(4,m) + backu(m)/BB(m)
            zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))
         enddo
         !$acc update device(zion(:,m0:mi))

         call n_T_u_group(m0,mi,flag_backf0=.true.)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))
         call ion_locator_at_larmor_points_partial(m0,mi) ! This may not be necessary, but just in case something goes wrong in the future

         allocate(target_norm_phase_vol(m0:mi))
         if(present(mesh_norm_dxdv)) then

            mesh_min_x = [upara_min, uperp_min, a_min, gq_min, gf_min]
            mesh_inv_dx = [1.0/dupara, 1.0/duperp, 1.0/da, 1.0/dgq, 1.0/dgf]
            mesh_dims = [dim_upara, dim_uperp, dim_a, dim_gq, dim_gf]

            !$acc parallel loop present(target_norm_phase_vol) private(ptl_data,mesh_idx)
            do m=m0,mi
               ptl_data(1) = (zion(4,m)*BB(m) -backu(m))/sqrt(backT(m)) ! upara
               ptl_data(2) = zion(6,m)*sqrt(2.0*BB(m)/backT(m)) ! uperp
               ptl_data(3) = zion(1,m) ! a
               ptl_data(4) = zion(2,m) ! gq
               ptl_data(5) = zion(3,m) ! gf

               !$acc loop seq
               do i=1,5
                  mesh_idx(i) = ceiling( (ptl_data(i)-mesh_min_x(i))*mesh_inv_dx(i) + 0.5 ) ! round operation to find nearst grid
               enddo

               target_norm_phase_vol(m) = 1.0-mesh_norm_dxdv(mesh_idx(1),mesh_idx(2),mesh_idx(3),mesh_idx(4),mesh_idx(5))
            enddo
         else
            !$acc parallel loop present(target_norm_phase_vol)
            do m=m0,mi
               target_norm_phase_vol(m) = 1.0
            enddo
         endif


         ! Set initial weights of ions
         ! zion(5,m) is "df"
         ! zion0(8,m) is "f_tot"
         ! zion0(6,m) is "phase volume"

         !$acc parallel loop present(zion,zion0,target_norm_phase_vol,back_f0_i,wi_init,wi_vpara,wi_vperp,wi_accel,wi_fdamp,UnNorm_theta_i,UnNorm_zeta_i)
         do m=m0,mi
            zion0(6,m) = target_norm_phase_vol(m)/gi_val(m) ! Phase space volume
            integral_f0 = zion0(6,m)*back_f0_i(m) !f_0
            target_w_tot = integral_f0

            zion0(8,m) = target_w_tot ! full-f
            zion(5,m) = zion0(8,m)-integral_f0 ! df
            zion(8,m) = integral_f0 ! f_0

            wi_init(m)  = zion(5,m)
            wi_vpara(m) = 0
            wi_vperp(m) = 0
            wi_accel(m) = 0
            wi_fdamp(m) = 0

            !$acc loop seq
            do i=1,5
               zion0(i,m)     = zion(i,m)
            enddo
            zion0(7,m)       = zion(7,m)
            UnNorm_theta_i(m)= zion(2,m)
            UnNorm_zeta_i(m) = zion(3,m)
         enddo
         acc_newer_D_zion=.true.

         deallocate(gi_val, rand_5d)
         deallocate(target_norm_phase_vol)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Generate_Ions_in_given_Phase_Space takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Generate_Ions_in_given_Phase_Space

      subroutine Generate_Electrons_in_given_Phase_Space(upara,uperp,a,gq,gf,this_N_newP, mesh_norm_dxdv,in_label_val)
         use global_parameters
         use particle_array
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         use particle_decomp
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         real(wp), intent(in), dimension(:) :: upara,uperp,a,gq,gf
         integer, intent(in) :: this_N_newP
         real(wp), intent(in), optional :: mesh_norm_dxdv(:,:,:,:,:) ! 5-D grid (upara,uperp,a,gq,gf) normalized dxdv by pre-occupying particles
         real(wp), intent(in), optional :: in_label_val

         real(wp) :: label_val
         integer m,m0,i, tot_N_newP
         real(wp) :: vpara, vperp, integral_f0, target_w_tot
         real(wp), allocatable :: ge_val(:)
         real(wp), allocatable :: rand_5d(:,:)
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max, upara_min, upara_max, uperp_min, uperp_max
         real(wp) :: Del_a, Del_gq, Del_gf, Del_upara, Del_uperp
         real(wp) :: da, dgq, dgf, dupara, duperp
         integer :: dim_a, dim_gq, dim_gf, dim_upara, dim_uperp
         integer :: minpts, ifail=1, ierror
         real(wp) :: ac, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian
         real(wp) time0, time1

         integer ij,ip
         real(wp) fip,r,dela
         real(wp) wz0,wz1,wp0,wt00,tw10,wp1,wt01,wt11,wt10
         real(wp) intp_phi, intp_FLA_Te0, intp_phi00
         real(wp), allocatable, dimension(:) :: AdiabTerm ! e*phi/T

         real(wp) :: mesh_min_x(5),mesh_inv_dx(5)
         real(wp) :: ptl_data(5)
         integer :: mesh_idx(5), mesh_dims(5)
         real(wp), allocatable :: target_norm_phase_vol(:)
         !$acc declare create(target_norm_phase_vol, AdiabTerm)

         call MPI_ALLREDUCE(this_N_newP, tot_N_newP,1,MPI_INTEGER, MPI_SUM, partd_comm,ierror)

         if(tot_N_newP<=0) then
            return
         endif

         label_val=0
         if(present(in_label_val)) label_val = in_label_val

#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         upara_min = minval(upara); upara_max = maxval(upara); dim_upara = size(upara)
         uperp_min = minval(uperp); uperp_max = maxval(uperp); dim_uperp = size(uperp)
         a_min = minval(a); a_max = maxval(a); dim_a = size(a)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         Del_upara = upara_max-upara_min;  dupara=Del_upara/(dim_upara-1)
         Del_uperp = uperp_max-uperp_min;  duperp=Del_uperp/(dim_uperp-1)
         Del_a = a_max-a_min;       da=Del_a/(dim_a-1)
         Del_gq = gq_max-gq_min;    dgq=Del_gq/(dim_gq-1)
         Del_gf = gf_max-gf_min;    dgf=Del_gf/(dim_gf-1)


         ! Generate
         m0=me+1
         me=me+this_N_newP

         allocate(ge_val(m0:me), rand_5d(5,m0:me))

         ! Prepare random numbers
         do m=m0,me
            do i=1,5
               rand_5d(i,m) = sprng()
            enddo
         enddo


         do m=m0,me
            zeon(1,m) = a_min + rand_5d(1,m)*Del_a
            zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,ac,dsda,-1,ifail)

            zeon(2,m) = modulo(gq_min + rand_5d(2,m)*Del_gq,2.0*pi)  ! gq
            zeon(3,m) = modulo(gf_min + rand_5d(3,m)*Del_gf,2.0*pi)  ! gf

            elabels(m) = label_val
         enddo
         !$acc update device(zeon(:,m0:me),elabels(m0:me))

         call n_T_u_group_e(m0,me)   ! get Te, note that back_f0_e cannot be calculated here because it requires zeon(4) and zeon(6)
         !$acc update self(backTe(m0:me))
         call esigetbf_partial(zeon,m0,me,zeondim)


         ! Uniform Loading
         do m=m0,me
            vpara = (upara_min + rand_5d(4,m)*Del_upara)*sqrt(backTe(m))
            vperp = (uperp_min + rand_5d(5,m)*Del_uperp)*sqrt(backTe(m))
            ge_val(m) = tot_N_newP/(Del_a*Del_gq*Del_gf*Del_upara*Del_uperp*backTe(m))/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)

            zeon(4,m) = msm0*vpara/(csc0*BB(m))
            zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))
         enddo
         !$acc update device(zeon(:,m0:me))

         call n_T_u_group_e(m0,me,flag_backf0=.true.)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))

         allocate(target_norm_phase_vol(m0:me))
         if(present(mesh_norm_dxdv)) then

            mesh_min_x = [upara_min, uperp_min, a_min, gq_min, gf_min]
            mesh_inv_dx = [1.0/dupara, 1.0/duperp, 1.0/da, 1.0/dgq, 1.0/dgf]
            mesh_dims = [dim_upara, dim_uperp, dim_a, dim_gq, dim_gf]

            !$acc parallel loop present(target_norm_phase_vol) private(ptl_data,mesh_idx)
            do m=m0,me
               ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
               ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))
               ptl_data(1) = zeon(4,m)*csc0*BB(m)*sqrt(m0ms/(msm0*backTe(m)))! upara
               ptl_data(2) = zeon(6,m)*sqrt(2.0*BB(m))/sqrt(msm0*backTe(m))  ! uperp
               ptl_data(3) = zeon(1,m) ! a
               ptl_data(4) = zeon(2,m) ! gq
               ptl_data(5) = zeon(3,m) ! gf

               !$acc loop seq
               do i=1,5
                  mesh_idx(i) = ceiling( (ptl_data(i)-mesh_min_x(i))*mesh_inv_dx(i) + 0.5 ) ! round operation to find nearst grid
               enddo

               target_norm_phase_vol(m) = 1.0-mesh_norm_dxdv(mesh_idx(1),mesh_idx(2),mesh_idx(3),mesh_idx(4),mesh_idx(5))
            enddo
         else
            !$acc parallel loop present(target_norm_phase_vol)
            do m=m0,me
               target_norm_phase_vol(m) = 1.0
            enddo
         endif


         ! Calculate Adiab term
         allocate(AdiabTerm(m0:me))
         call electron_locator_partial(m0,me) ! This is required to calculate wpeon, wteon0 and etc
         !$acc parallel loop present(AdiabTerm, wzeon,jteon0,wpeon,wteon0,jteon1,wteon1)
         do m=m0,me
            wz1=wzeon(m)
            wz0=1.0_wp-wz1

            ij=jteon0(m)
            wp0=1.0_wp-wpeon(m)
            wt00=1.0_wp-wteon0(m)
            intp_phi = wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0 = wp0*wt00*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=ij+1
            wt10=1.0_wp-wt00
            intp_phi=intp_phi+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0=intp_FLA_Te0+wp0*wt10*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=jteon1(m)
            wp1=1.0_wp-wp0
            wt01=1.0_wp-wteon1(m)
            intp_phi=intp_phi+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0=intp_FLA_Te0+wp1*wt01*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=ij+1
            wt11=1.0_wp-wt01
            intp_phi=intp_phi+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons
            intp_FLA_Te0=intp_FLA_Te0+wp1*wt11*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

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


         ! Set initial weights of ions
         ! zeon(5,m) is "df"
         ! zeon(8,m) is "f_tot"
         ! zeon0(6,m) is "phase volume"

         !$acc parallel loop present(zeon,zeon0,target_norm_phase_vol,back_f0_e,we_init,we_vpara,we_vperp,we_accel,we_fdamp,UnNorm_theta_e,UnNorm_zeta_e)
         do m=m0,me
            ! zeon(5,m) is "df" or "dh"
            ! zeon(8,m) is "full-f or full-h"
            ! zeon(9,m) is "f_AD"
            ! zeon0(9,m) is "f_tot"
            ! zeon0(6,m) is "phase volume"

            zeon0(6,m) = target_norm_phase_vol(m)/ge_val(m) ! Phase space volume
            integral_f0 = zeon0(6,m)*back_f0_e(m) !f_0

            target_w_tot = integral_f0

            ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
            zeon(9,m) = integral_f0*AdiabTerm(m) ! f_AD

            zeon(8,m) = target_w_tot -int(ifelectron/2)*zeon(9,m)! full-f (ifelectron<2) or full-h (ifelectron==2)
            zeon(5,m)  = zeon(8,m)-integral_f0 ! df or dh
            zeon0(9,m) = target_w_tot ! f_tot

            we_init(m)  = zeon(5,m)
            we_vpara(m) = 0
            we_vperp(m) = 0
            we_accel(m) = 0
            we_adiab(m) = 0
            we_fdamp(m) = 0

            !$acc loop seq
            do i=1,5
               zeon0(i,m)     = zeon(i,m)
            enddo
            zeon0(7,m)       = zeon(7,m)
            zeon0(8,m)       = zeon(8,m)
            UnNorm_theta_e(m)= zeon(2,m)
            UnNorm_zeta_e(m) = zeon(3,m)
         enddo
         acc_newer_D_zeon=.true.

         deallocate(ge_val, rand_5d)
         deallocate(target_norm_phase_vol)
         deallocate(AdiabTerm)


#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Generate_Electrons_in_given_Phase_Space takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Generate_Electrons_in_given_Phase_Space


      subroutine Generate_Electrons_in_given_Phase_Space_PolarVelCoord(ur,uq,a,gq,gf,this_N_newP, mesh_norm_dxdv,in_label_val)
         use global_parameters
         use particle_array
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         use particle_decomp
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         real(wp), intent(in), dimension(:) :: ur,uq,a,gq,gf
         integer, intent(in) :: this_N_newP
         real(wp), intent(in), optional :: mesh_norm_dxdv(:,:,:,:,:) ! 5-D grid (ur,uq,a,gq,gf) normalized dxdv by pre-occupying particles
         real(wp), intent(in), optional :: in_label_val

         real(wp) :: label_val, upara, uperp, tmp_ur, tmp_uq
         integer m,m0,i, tot_N_newP
         real(wp) :: vpara, vperp, integral_f0, target_w_tot
         real(wp), allocatable :: ge_val(:)
         real(wp), allocatable :: rand_5d(:,:)
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max, ur_min, ur_max, uq_min, uq_max
         real(wp) :: Del_a, Del_gq, Del_gf, Del_ur, Del_uq
         real(wp) :: da, dgq, dgf, dur, duq
         integer :: dim_a, dim_gq, dim_gf, dim_ur, dim_uq
         integer :: minpts, ifail=1, ierror
         real(wp) :: ac, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian
         real(wp) time0, time1

         integer ij,ip
         real(wp) fip,r,dela
         real(wp) wz0,wz1,wp0,wt00,tw10,wp1,wt01,wt11,wt10
         real(wp) intp_phi, intp_FLA_Te0, intp_phi00
         real(wp), allocatable, dimension(:) :: AdiabTerm ! e*phi/T

         real(wp) :: mesh_min_x(5),mesh_inv_dx(5)
         real(wp) :: ptl_data(5)
         integer :: mesh_idx(5), mesh_dims(5)
         real(wp), allocatable :: target_norm_phase_vol(:)
         !$acc declare create(target_norm_phase_vol, AdiabTerm)

         call MPI_ALLREDUCE(this_N_newP, tot_N_newP,1,MPI_INTEGER, MPI_SUM, partd_comm,ierror)

         if(tot_N_newP<=0) then
            return
         endif

         label_val=0
         if(present(in_label_val)) label_val = in_label_val

#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

         ur_min = minval(ur); ur_max = maxval(ur); dim_ur = size(ur)
         uq_min = minval(uq); uq_max = maxval(uq); dim_uq = size(uq)
         a_min = minval(a); a_max = maxval(a); dim_a = size(a)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         Del_ur = ur_max-ur_min;  dur=Del_ur/(dim_ur-1)
         Del_uq = uq_max-uq_min;  duq=Del_uq/(dim_uq-1)
         Del_a = a_max-a_min;       da=Del_a/(dim_a-1)
         Del_gq = gq_max-gq_min;    dgq=Del_gq/(dim_gq-1)
         Del_gf = gf_max-gf_min;    dgf=Del_gf/(dim_gf-1)


         ! Generate
         m0=me+1
         me=me+this_N_newP

         allocate(ge_val(m0:me), rand_5d(5,m0:me))

         ! Prepare random numbers
         do m=m0,me
            do i=1,5
               rand_5d(i,m) = sprng()
            enddo
         enddo


         do m=m0,me
            zeon(1,m) = a_min + rand_5d(1,m)*Del_a
            zeon(7,m)=D01AHF(0.0_wp,zeon(1,m),error_8,minpts,ac,dsda,-1,ifail)

            zeon(2,m) = modulo(gq_min + rand_5d(2,m)*Del_gq,2.0*pi)  ! gq
            zeon(3,m) = modulo(gf_min + rand_5d(3,m)*Del_gf,2.0*pi)  ! gf

            elabels(m) = label_val
         enddo
         !$acc update device(zeon(:,m0:me),elabels(m0:me))

         call n_T_u_group_e(m0,me)   ! get Te, note that back_f0_e cannot be calculated here because it requires zeon(4) and zeon(6)
         !$acc update self(backTe(m0:me))
         call esigetbf_partial(zeon,m0,me,zeondim)


         ! Uniform Loading
         do m=m0,me
            ! Find sampling points uniformly distributed within an annulus in the poloar coordinates (vpara,vperp)
            ! r=sqrt(a1^2+rand*(a2^2-a1^2));
            tmp_ur = sqrt(ur_min**2.0 + rand_5d(4,m)*(ur_max**2.0-ur_min**2.0))
            tmp_uq = (uq_min + rand_5d(5,m)*Del_uq)

            vpara = (tmp_ur*cos(tmp_uq))*sqrt(backTe(m))
            vperp = (tmp_ur*sin(tmp_uq))*sqrt(backTe(m))

            ! Note that ge_val (marker distribution) should include the Jacobian from the polar coordinate (here, tmp_ur)
            ge_val(m) = tot_N_newP/(Del_a*Del_gq*Del_gf*(0.5*pi*(ur_max**2.0-ur_min**2.0))*backTe(m))/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)

            zeon(4,m) = msm0*vpara/(csc0*BB(m))
            zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))
         enddo
         !$acc update device(zeon(:,m0:me))

         call n_T_u_group_e(m0,me,flag_backf0=.true.)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))

         allocate(target_norm_phase_vol(m0:me))
         if(present(mesh_norm_dxdv)) then

            mesh_min_x = [ur_min, uq_min, a_min, gq_min, gf_min]
            mesh_inv_dx = [1.0/dur, 1.0/duq, 1.0/da, 1.0/dgq, 1.0/dgf]
            mesh_dims = [dim_ur, dim_uq, dim_a, dim_gq, dim_gf]

            !$acc parallel loop present(target_norm_phase_vol) private(ptl_data,mesh_idx)
            do m=m0,me
               ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
               ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))
               upara = zeon(4,m)*csc0*BB(m)*sqrt(m0ms/(msm0*backTe(m)))! upara
               uperp = zeon(6,m)*sqrt(2.0*BB(m))/sqrt(msm0*backTe(m))  ! uperp

               ptl_data(1) = sqrt(upara**2.0+uperp**2.0) ! ur (raidus in the polar coordinate)
               ptl_data(2) = acos(upara/ptl_data(1)) ! uq (angle in the polar coordinate)
               ptl_data(3) = zeon(1,m) ! a
               ptl_data(4) = zeon(2,m) ! gq
               ptl_data(5) = zeon(3,m) ! gf

               !$acc loop seq
               do i=1,5
                  mesh_idx(i) = ceiling( (ptl_data(i)-mesh_min_x(i))*mesh_inv_dx(i) + 0.5 ) ! round operation to find nearst grid
               enddo

               target_norm_phase_vol(m) = 1.0-mesh_norm_dxdv(mesh_idx(1),mesh_idx(2),mesh_idx(3),mesh_idx(4),mesh_idx(5))
            enddo
         else
            !$acc parallel loop present(target_norm_phase_vol)
            do m=m0,me
               target_norm_phase_vol(m) = 1.0
            enddo
         endif


         ! Calculate Adiab term
         allocate(AdiabTerm(m0:me))
         call electron_locator_partial(m0,me) ! This is required to calculate wpeon, wteon0 and etc
         !$acc parallel loop present(AdiabTerm, wzeon,jteon0,wpeon,wteon0,jteon1,wteon1)
         do m=m0,me
            wz1=wzeon(m)
            wz0=1.0_wp-wz1

            ij=jteon0(m)
            wp0=1.0_wp-wpeon(m)
            wt00=1.0_wp-wteon0(m)
            intp_phi = wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0 = wp0*wt00*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=ij+1
            wt10=1.0_wp-wt00
            intp_phi=intp_phi+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0=intp_FLA_Te0+wp0*wt10*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=jteon1(m)
            wp1=1.0_wp-wp0
            wt01=1.0_wp-wteon1(m)
            intp_phi=intp_phi+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))
            intp_FLA_Te0=intp_FLA_Te0+wp1*wt01*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

            ij=ij+1
            wt11=1.0_wp-wt01
            intp_phi=intp_phi+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons
            intp_FLA_Te0=intp_FLA_Te0+wp1*wt11*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

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


         ! Set initial weights of ions
         ! zeon(5,m) is "df"
         ! zeon(8,m) is "f_tot"
         ! zeon0(6,m) is "phase volume"

         !$acc parallel loop present(zeon,zeon0,target_norm_phase_vol,back_f0_e,we_init,we_vpara,we_vperp,we_accel,we_fdamp,UnNorm_theta_e,UnNorm_zeta_e)
         do m=m0,me
            ! zeon(5,m) is "df" or "dh"
            ! zeon(8,m) is "full-f or full-h"
            ! zeon(9,m) is "f_AD"
            ! zeon0(9,m) is "f_tot"
            ! zeon0(6,m) is "phase volume"

            zeon0(6,m) = target_norm_phase_vol(m)/ge_val(m) ! Phase space volume
            integral_f0 = zeon0(6,m)*back_f0_e(m) !f_0

            target_w_tot = integral_f0

            ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
            zeon(9,m) = integral_f0*AdiabTerm(m) ! f_AD

            zeon(8,m) = target_w_tot -int(ifelectron/2)*zeon(9,m)! full-f (ifelectron<2) or full-h (ifelectron==2)
            zeon(5,m)  = zeon(8,m)-integral_f0 ! df or dh
            zeon0(9,m) = target_w_tot ! f_tot

            we_init(m)  = zeon(5,m)
            we_vpara(m) = 0
            we_vperp(m) = 0
            we_accel(m) = 0
            we_adiab(m) = 0
            we_fdamp(m) = 0

            !$acc loop seq
            do i=1,5
               zeon0(i,m)     = zeon(i,m)
            enddo
            zeon0(7,m)       = zeon(7,m)
            zeon0(8,m)       = zeon(8,m)
            UnNorm_theta_e(m)= zeon(2,m)
            UnNorm_zeta_e(m) = zeon(3,m)
         enddo
         acc_newer_D_zeon=.true.

         deallocate(ge_val, rand_5d)
         deallocate(target_norm_phase_vol)
         deallocate(AdiabTerm)


#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Generate_Electrons_in_given_Phase_Space_PolarVelCoord takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Generate_Electrons_in_given_Phase_Space_PolarVelCoord

      subroutine Generate_Ions_in_given_Phase_Space_PolarVelCoord_AdapTgtUr(a,gq,gf,tgt_ur_bound,mesh_occp_MarkerDen,mesh5D_norm_occp_dxdv,in_label_val)
         use global_parameters
         use particle_array
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         use particle_decomp
         use Structured_Grid
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         real(wp), intent(in), dimension(:) :: a,gq,gf
         real(wp), intent(in) :: tgt_ur_bound(:,:,:,:) ! (dim_a,dim_gq,dim_gf,2); upper ur and lower ur bound on 3-D grid (a,gq,gf)
         real(wp), intent(in), optional, dimension(:,:,:,:,:) :: mesh_occp_MarkerDen, mesh5D_norm_occp_dxdv ! 5-D grid (sr,sq,a,gq,gf) normalized MarkerDensity and dxdv (phase volume) occupied by present particles
         real(wp), intent(in), optional :: in_label_val
         ! where sr (speed) refers normalized ur regarding tgt_ur_bound:
         ! s(sr,sq,a,gq,gf) = (u(sr,sq,a,gq,gf)-ur_LowerBound(a,gq,gf))/(ur_UpperBound(a,gq,gf)-ur_LowerBound(a,gq,gf))
         ! sr(radius) and sq(angle) are polar coordinates for velocity space
         ! Particles will be generated in the domain of s=[0,1] on every spatial node
         integer :: dim_sr, dim_sq, srid, sqid
         real(wp) :: sr_min, sr_max, sq_min, sq_max
         real(wp) :: tgt_cell_Nm, sr0,sr1,ur0,ur1,sq0,sq1
         real(wp),allocatable :: sr_1D(:), sq_1D(:)
         real(wp) :: Del_sr, Del_sq
         real(wp) ini_MarkerDen_mype, ini_MarkerDen_npartd
         integer :: BLK_SZa, BLK_SZq  ! Block size (number of indices) into a, gq direction
         integer :: aid,qid,fid, loc_max_aid, loc_max_qid, loc_Nsamp

         real(wp) :: label_val, upara, uperp, tmp_ur, tmp_uq
         integer m,m0,i, tot_N_newP
         real(wp) :: gi_val, vpara, vperp, integral_f0, target_w_tot
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max
         real(wp) :: Del_a, Del_gq, Del_gf
         real(wp) :: da, dgq, dgf, dsr, dsq
         integer :: dim_a, dim_gq, dim_gf
         integer :: minpts, ifail=1, ierror
         real(wp) :: ac, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian
         real(wp) time0, time1

         real(wp) loc_max_ur, loc_min_ur, loc_Del_a, loc_Del_q, loc_Del_f, loc_gf_min, loc_gf_max
         real(wp) ptl_a, ptl_gq, ptl_sr,ptl_sq, ptl_ur, inv_da, inv_dgq,inv_dsr,inv_dsq, ws0,ws1
         integer tmp_aid,k,ptl_aid,ptl_gqid, ptl_srid, ptl_sqid, Bid
         real(wp) P_safety_q
         real(wp) ptl_ur_bound(2)
         integer, allocatable :: num_BLKs,BLKs_aid(:),BLKs_gqid(:),BLKs_gfid(:), BLKs_m0(:), BLKs_m1(:)
         integer, allocatable :: BLKs_Nsamp_mype(:), BLKs_Nsamp_npartd(:),BLKs_NnewP_mype(:), BLKs_NnewP_npartd(:)
         real(wp), allocatable :: BLKs_d3x(:), BLKs_d2u(:) ! spatial area and velocity area of each block
         real(wp) occp_MarkerDen,amp_factor, mean_occp_MarkerDen, mean_norm_occp_dxdv
         character(len=128) filename

         label_val=0
         if(present(in_label_val)) label_val = in_label_val


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif


         sr_min = 0.0; sr_max = 1.0;
         sq_min = 0.0; sq_max = pi;

         if(present(mesh_occp_MarkerDen)) then
            dim_sr = size(mesh_occp_MarkerDen,1)
            dim_sq = size(mesh_occp_MarkerDen,2)
            if(present(mesh5D_norm_occp_dxdv)) then
               if(dim_sr/=size(mesh5D_norm_occp_dxdv,1)) then
                  write(stdout,*) "Error!! mesh5D_MarkerDen's dim_sr=",dim_sr," is different from mesh5D_dxdv's dim_s=",size(mesh5D_norm_occp_dxdv,1)
                  call flush(stdout)
                  stop
               endif
               if(dim_sq/=size(mesh5D_norm_occp_dxdv,2)) then
                  write(stdout,*) "Error!! mesh5D_MarkerDen's dim_sq=",dim_sq," is different from mesh5D_dxdv's dim_q=",size(mesh5D_norm_occp_dxdv,2)
                  call flush(stdout)
                  stop
               endif
            endif
         elseif(present(mesh5D_norm_occp_dxdv)) then
            dim_sr = size(mesh5D_norm_occp_dxdv,1)
            dim_sq = size(mesh5D_norm_occp_dxdv,2)
         endif


         a_min = minval(a); a_max = maxval(a); dim_a = size(a)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         Del_sr = sr_max-sr_min;    dsr=Del_sr/(dim_sr-1)
         Del_sq = sq_max-sq_min;    dsq=Del_sq/(dim_sq-1)
         Del_a = a_max-a_min;       da=Del_a/(dim_a-1)
         Del_gq = gq_max-gq_min;    dgq=Del_gq/(dim_gq-1)
         Del_gf = gf_max-gf_min;    dgf=Del_gf/(dim_gf-1)


         allocate(sr_1D(dim_sr),sq_1D(dim_sq))

         do srid=1,dim_sr
            sr_1D(srid)=sr_min+(srid-1)*dsr
         enddo
         do sqid=1,dim_sq
            sq_1D(sqid)=sq_min+(sqid-1)*dsq
         enddo

         ini_MarkerDen_mype = ini_me/((a1-a0)*(2*pi)*(zetamax-zetamin)*0.5*pi*umax**2.0)
         ini_MarkerDen_npartd = ini_me_npartd/((a1-a0)*(2*pi)*(zetamax-zetamin)*0.5*pi*umax**2.0)



         inv_da = 1.0/da; inv_dgq=1.0/dgq; inv_dsr=1.0/dsr; inv_dsq=1.0/dsq

         ! BLK_SZa = PtlSrcCV%dim_a-1
         ! BLK_SZq = PtlSrcCV%dim_gq-1
         BLK_SZa = 7
         BLK_SZq = 20

         num_BLKs = dim_gf*ceiling(real(dim_a-1)/BLK_SZa)*ceiling(real(dim_gq-1)/BLK_SZq)
         allocate(BLKs_aid(num_BLKs),BLKs_gqid(num_BLKs),BLKs_gfid(num_BLKs),BLKs_m0(num_BLKs), BLKs_m1(num_BLKs))
         allocate(BLKs_Nsamp_mype(num_BLKs), BLKs_Nsamp_npartd(num_BLKs),BLKs_NnewP_mype(num_BLKs),BLKs_NnewP_npartd(num_BLKs))
         allocate(BLKs_d3x(num_BLKs),BLKs_d2u(num_BLKs))


         Bid=0
         do fid=1,dim_gf
            do aid=1,dim_a-1,BLK_SZa
               do qid=1,dim_gq-1,BLK_SZq
                  Bid = Bid + 1
                  BLKs_aid(Bid)=aid
                  BLKs_gqid(Bid)=qid
                  BLKs_gfid(Bid)=fid
               enddo
            enddo
         enddo

         if(Bid/=num_BLKs) then
            write(stdout,'((A,I0,A,I0))') "Error!! Something is Wrong!! After finding block idices, Bid=", Bid, "; num_BLKs=", num_BLKs
            stop
         endif

         m0 = mi+1

         do Bid=1,num_BLKs

            BLKs_m0(Bid) = mi+1

            aid = BLKs_aid(Bid)
            qid = BLKs_gqid(Bid)
            fid = BLKs_gfid(Bid)

            loc_max_aid = min(dim_a,aid+BLK_SZa)
            loc_Del_a = a(loc_max_aid)-a(aid)
            loc_max_qid = min(dim_gq,qid+BLK_SZq)
            loc_Del_q = gq(loc_max_qid)-gq(qid)

            loc_gf_max = min(gf(dim_gf),gf(fid)+0.5*dgf)
            loc_gf_min = max(gf(1),gf(fid)-0.5*dgf)
            loc_Del_f = loc_gf_max - loc_gf_min

            BLKs_d3x(Bid) = loc_Del_f*loc_Del_a*loc_Del_q

            loc_min_ur = minval(tgt_ur_bound(aid:loc_max_aid, qid:loc_max_qid, fid,1)) ! lower bound
            loc_max_ur = maxval(tgt_ur_bound(aid:loc_max_aid, qid:loc_max_qid, fid,2)) ! upper bound

            BLKs_d2u(Bid) = 0.5*pi*(loc_max_ur**2.0-loc_min_ur**2.0)

            loc_Nsamp = ceiling(ini_MarkerDen_mype*BLKs_d3x(Bid)*BLKs_d2u(Bid) +sprng())
            BLKs_Nsamp_mype(Bid) = loc_Nsamp


            BLKs_NnewP_mype(Bid) = 0
            do k=1,loc_Nsamp
               ptl_a  = a(aid) +loc_Del_a*sprng()
               ptl_gq = gq(qid)+loc_Del_q*sprng() ! NOTE!! This is gq value at either zeta=zetamin or zeta=zetamax

               ptl_aid  = ceiling((ptl_a-a_min)*inv_da+0.5) ! Find nearst grid
               ptl_gqid = ceiling((ptl_gq-gq_min)*inv_dgq+0.5) ! Find nearst grid

               ! ur lower(1) and upper(2) bound at ptl position
               ptl_ur_bound = tgt_ur_bound(ptl_aid,ptl_gqid,fid,:)

               ! Find sampling points uniformly distributed within an annulus in the poloar coordinates (vpara,vperp)
               ! r=sqrt(a1^2+rand*(a2^2-a1^2));
               ptl_ur = sqrt(loc_min_ur**2.0 + sprng()*(loc_max_ur**2.0-loc_min_ur**2.0))


               if( (ptl_ur-ptl_ur_bound(1))*(ptl_ur-ptl_ur_bound(2)) < 0 ) then
                  ! sr (speed) refers normalized ur regarding tgt_ur_bound:
                  ! sr(sr,sq,a,gq,gf) = (ur(sr,sq,a,gq,gf)-ur_LowerBound(a,gq,gf))/(ur_UpperBound(a,gq,gf)-ur_LowerBound(a,gq,gf))
                  ptl_sr = (ptl_ur-ptl_ur_bound(1))/(ptl_ur_bound(2)-ptl_ur_bound(1))
                  ptl_srid = ceiling((ptl_sr-sr_min)*inv_dsr+0.5)

                  ptl_sq = sq_min + sprng()*(sq_max-sq_min)
                  ptl_sqid = ceiling((ptl_sq-sq_min)*inv_dsq+0.5)

                  !accept this particle
                  BLKs_NnewP_mype(Bid) = BLKs_NnewP_mype(Bid) +1
                  mi = mi+1

                  ! This is only a part of zion0(6,m); it will be recalulated by multiplying 1/ge_val later
                  if(present(mesh5D_norm_occp_dxdv)) then
                     zion0(6,mi) = 1.0 -mesh5D_norm_occp_dxdv(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid)
                  else
                     zion0(6,mi) = 1.0
                  endif

                  if(present(mesh_occp_MarkerDen)) then

                     sq0 = max(sq_min, sq_1D(ptl_sqid)-0.5*dsq)
                     sq1 = min(sq_max, sq_1D(ptl_sqid)+0.5*dsq)

                     sr0 = max(sr_min, sr_1D(ptl_srid)-0.5*dsr)
                     sr1 = min(sr_max, sr_1D(ptl_srid)+0.5*dsr)
                     ur0 = ptl_ur_bound(1) + sr0*(ptl_ur_bound(2)-ptl_ur_bound(1))
                     ur1 = ptl_ur_bound(1) + sr1*(ptl_ur_bound(2)-ptl_ur_bound(1))

                     ! ! mesh_cell_area =  Cell_dVol_x3*0.5*(ur1**2.0-ur0**2.0)*(sq1-sq0)
                     tgt_cell_Nm = ini_MarkerDen_mype*(min(a1,a(ptl_aid)+0.5*da)-max(a0,a(ptl_aid)-0.5*da))*dgq*(dgf)*0.5*(ur1**2.0-ur0**2.0)*(sq1-sq0)

                     ! occp_MarkerDen = (1.0/npartdom)*mesh_occp_MarkerDen(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid)
                     occp_MarkerDen = mesh_occp_MarkerDen(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid)

                     if(occp_MarkerDen>0.95*ini_MarkerDen_npartd .and. occp_MarkerDen<=ini_MarkerDen_npartd) then
                        occp_MarkerDen=0.95*ini_MarkerDen_npartd ! To prevent too close occp_MarkerDen~ini_MarkerDen which could be troublesome
                     endif

                     if(sprng()*ini_MarkerDen_npartd <occp_MarkerDen) then
                        BLKs_NnewP_mype(Bid) = BLKs_NnewP_mype(Bid) -1
                        mi=mi-1
                        cycle ! reject this particle
                     endif
                     ! Note that denominator is always smaller than 1 otherwise it should be rejected by the previous checking routine
                     zion0(6,mi) = zion0(6,mi)/(1.0-occp_MarkerDen/ini_MarkerDen_npartd)
                  endif


                  zion(1,mi) = ptl_a
                  zion(7,mi) = D01AHF(0.0_wp,zion(1,mi),error_8,minpts,ac,dsda,-1,ifail)

                  ! Interpolate safety factor
                  tmp_aid = max(1,min(int((zion(1,mi)-PSCV%a_min)*PSCV%inv_da) +1,PSCV%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
                  ws1 = (zion(1,mi)-PSCV%a_1D(tmp_aid))*PSCV%inv_da
                  ws0 = (1.0-ws1)
                  P_safety_q = ws0*PSCV%safety_q_1D(tmp_aid) + ws1*PSCV%safety_q_1D(tmp_aid+1)

                  zion(3,mi) = modulo(loc_gf_min + loc_Del_f*sprng(), 2.0*pi)
                  ! zion(2,mi) = modulo(ptl_gq + (gf(fid)-zion(3,mi))/P_safety_q, 2.0*pi)
                  zion(2,mi) = modulo(ptl_gq + (zion(3,mi)-gf(fid))/P_safety_q, 2.0*pi)

                  upara = (ptl_ur*cos(ptl_sq))
                  uperp = (ptl_ur*sin(ptl_sq))

                  ! Note that zion(4) and zion(6) are temporary used to carray upara and uperp values
                  ! The real zion values will be calculated later
                  zion(4,mi) = upara
                  zion(6,mi) = uperp

                  ilabels(mi) = label_val
               endif
            enddo ! End k=1,loc_Nsamp

            BLKs_m1(Bid) = mi

         enddo ! End BLK loop
         !$acc update device(zion(:,m0:mi),ilabels(m0:mi))
         call MPI_ALLREDUCE(BLKs_NnewP_mype, BLKs_NnewP_npartd, num_BLKs, MPI_INTEGER, MPI_SUM, partd_comm,ierror)
         call MPI_ALLREDUCE(BLKs_Nsamp_mype, BLKs_Nsamp_npartd, num_BLKs, MPI_INTEGER, MPI_SUM, partd_comm,ierror)

         if(mype==0) then
            write(stdout,'(3(A,I0),A,F8.3)') "num_BLKs=",num_BLKs, "; tot_BLKs_NnewP_npartd=",sum(BLKs_NnewP_npartd),"; tot_BLKs_Nsamp_npartd=",sum(BLKs_Nsamp_npartd),"; Acceptance Rate=",real(sum(BLKs_NnewP_npartd))/sum(BLKs_Nsamp_npartd)
            write(stdout,'(3(A,I0))') "tot_BLKs_NnewP_mype=",sum(BLKs_NnewP_mype),"; mi=",mi,"; mimax=",mimax
            call flush(stdout)
         endif

         call n_T_u_group(m0,mi)   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(backT(m0:mi),backu(m0:mi))
         call esigetbf_partial(zion,m0,mi,ziondim)

         do Bid=1,num_BLKs
            do m=BLKs_m0(Bid),BLKs_m1(Bid)
               upara = zion(4,m) ! upara value is carried by zion(4) from the previous do loop
               uperp = zion(6,m) ! uperp value is carried by zion(6) from the previous do loop

               vpara = upara*sqrt(backT(m))
               vperp = uperp*sqrt(backT(m))

               zion(4,m) = (vpara+backu(m))/BB(m)
               zion(6,m) = vperp*sqrt(1.0/(2.0*BB(m)))

               ! Note that ge_val (marker distribution) should include the Jacobian from the polar coordinate (here, vperp)
               ! ge_val(m) = tot_N_newP/(Del_a*Del_gq*Del_gf*(0.5*pi*(ur_max**2.0-ur_min**2.0))*backT(m))/(Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)
               ! ge_val = BLKs_NnewP_npartd(Bid)/(BLKs_d3x(Bid)*BLKs_d2u(Bid)*backT(m)*Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)
               gi_val = BLKs_Nsamp_npartd(Bid)/(BLKs_d3x(Bid)*BLKs_d2u(Bid)*backT(m)*Jacobian(zion(1,m),zion(2,m))*2.0*pi*vperp)
               zion0(6,m) = zion0(6,m)/gi_val ! Phase space volume
            enddo
         enddo
         !$acc update device(zion(:,m0:mi))
         !$acc update device(zion0(6,m0:mi))
         call n_T_u_group(m0,mi,flag_backf0=.true.)   ! This is to get back_f0_i inside. (which requires zion(4) and zion(6))


         ! Set initial weights of ions
         ! zion(5,m) is "df"
         ! zion0(8,m) is "f_tot"
         ! zion0(6,m) is "phase volume"

         !$acc parallel loop present(zion,zion0,back_f0_i,wi_init,wi_vpara,wi_vperp,wi_accel,wi_fdamp,UnNorm_theta_i,UnNorm_zeta_i)
         do m=m0,mi
            ! zion(5,m) is "df"
            ! zion0(8,m) is "f_tot"
            ! zion(8,m) is "f_0"
            ! zion0(6,m) is "phase volume"

            integral_f0 = zion0(6,m)*back_f0_i(m) !f_0
            zion0(8,m) = integral_f0 ! full-f
            zion(5,m) = 0.0! df
            zion(8,m) = integral_f0 ! f_0

            wi_init(m)  = zion(5,m)
            wi_vpara(m) = 0
            wi_vperp(m) = 0
            wi_accel(m) = 0
            wi_fdamp(m) = 0

            !$acc loop seq
            do i=1,5
               zion0(i,m)     = zion(i,m)
            enddo
            zion0(7,m)       = zion(7,m)
            UnNorm_theta_i(m)= zion(2,m)
            UnNorm_zeta_i(m) = zion(3,m)
         enddo
         acc_newer_D_zion=.true.



         deallocate(BLKs_aid,BLKs_gqid,BLKs_gfid,BLKs_m0,BLKs_m1)
         deallocate(BLKs_Nsamp_mype, BLKs_Nsamp_npartd,BLKs_NnewP_mype,BLKs_NnewP_npartd)
         deallocate(BLKs_d3x,BLKs_d2u)
         deallocate(sr_1D,sq_1D)


#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Generate_Ions_in_given_Phase_Space_PolarVelCoord_AdapTgtUr takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Generate_Ions_in_given_Phase_Space_PolarVelCoord_AdapTgtUr


      subroutine Generate_Electrons_in_given_Phase_Space_PolarVelCoord_AdapTgtUr(a,gq,gf,tgt_ur_bound,mesh_occp_MarkerDen,mesh5D_norm_occp_dxdv,in_label_val)
         use global_parameters
         use particle_array
         use field_array
         use esi_interface
         use mpi_defs
         use function_gg_module
         use particle_decomp
         use Structured_Grid
         implicit none

#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng_f.h"

         real(wp), intent(in), dimension(:) :: a,gq,gf
         real(wp), intent(in) :: tgt_ur_bound(:,:,:,:) ! (dim_a,dim_gq,dim_gf,2); upper ur and lower ur bound on 3-D grid (a,gq,gf)
         real(wp), intent(in), optional, dimension(:,:,:,:,:) :: mesh_occp_MarkerDen, mesh5D_norm_occp_dxdv ! 5-D grid (sr,sq,a,gq,gf) normalized MarkerDensity and dxdv (phase volume) occupied by present particles
         real(wp), intent(in), optional :: in_label_val
         ! where sr (speed) refers normalized ur regarding tgt_ur_bound:
         ! s(sr,sq,a,gq,gf) = (u(sr,sq,a,gq,gf)-ur_LowerBound(a,gq,gf))/(ur_UpperBound(a,gq,gf)-ur_LowerBound(a,gq,gf))
         ! sr(radius) and sq(angle) are polar coordinates for velocity space
         ! Particles will be generated in the domain of s=[0,1] on every spatial node
         integer :: dim_sr, dim_sq, srid, sqid
         real(wp) :: sr_min, sr_max, sq_min, sq_max
         real(wp) :: tgt_cell_Nm, sr0,sr1,ur0,ur1,sq0,sq1
         real(wp),allocatable :: sr_1D(:), sq_1D(:)
         real(wp) :: Del_sr, Del_sq
         real(wp) ini_MarkerDen_mype, ini_MarkerDen_npartd
         integer :: BLK_SZa, BLK_SZq  ! Block size (number of indices) into a, gq direction
         integer :: aid,qid,fid, loc_max_aid, loc_max_qid, loc_Nsamp

         real(wp) :: label_val, upara, uperp, tmp_ur, tmp_uq
         integer m,m0,i, tot_N_newP
         real(wp) :: ge_val, vpara, vperp, integral_f0, target_w_tot
         real(wp) :: a_min, a_max, gq_min, gq_max, gf_min, gf_max
         real(wp) :: Del_a, Del_gq, Del_gf
         real(wp) :: da, dgq, dgf, dsr, dsq
         integer :: dim_a, dim_gq, dim_gf
         integer :: minpts, ifail=1, ierror
         real(wp) :: ac, error_8=1e-5
         real(wp), external:: dsda, D01AHF,jacobian
         real(wp) time0, time1

         integer ij,ip
         real(wp) fip,r,dela
         real(wp) wz0,wz1,wp0,wt00,tw10,wp1,wt01,wt11,wt10
         real(wp) intp_phi, intp_FLA_Te0, intp_phi00
         real(wp) Te_inverse, intp_phi_V, AD
         ! real(wp), allocatable, dimension(:) :: AdiabTerm ! e*phi/T

         real(wp) loc_max_ur, loc_min_ur, loc_Del_a, loc_Del_q, loc_Del_f, loc_gf_min, loc_gf_max
         real(wp) ptl_a, ptl_gq, ptl_sr,ptl_sq, ptl_ur, inv_da, inv_dgq,inv_dsr,inv_dsq, ws0,ws1
         integer tmp_aid,k,ptl_aid,ptl_gqid, ptl_srid, ptl_sqid, Bid
         real(wp) P_safety_q
         real(wp) ptl_ur_bound(2)
         integer, allocatable :: num_BLKs,BLKs_aid(:),BLKs_gqid(:),BLKs_gfid(:), BLKs_m0(:), BLKs_m1(:)
         integer, allocatable :: BLKs_Nsamp_mype(:), BLKs_Nsamp_npartd(:),BLKs_NnewP_mype(:), BLKs_NnewP_npartd(:)
         real(wp), allocatable :: BLKs_d3x(:), BLKs_d2u(:) ! spatial area and velocity area of each block
         real(wp) occp_MarkerDen,amp_factor, mean_occp_MarkerDen, mean_norm_occp_dxdv
         character(len=128) filename

         !!! $acc declare create(AdiabTerm)



         label_val=0
         if(present(in_label_val)) label_val = in_label_val


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif



         ! sr_min = 0.0; sr_max = 1.0; dim_sr = 6
         ! sq_min = 0.0; sq_max = pi;  dim_sq = 10
         sr_min = 0.0; sr_max = 1.0;
         sq_min = 0.0; sq_max = pi;

         if(present(mesh_occp_MarkerDen)) then
            dim_sr = size(mesh_occp_MarkerDen,1)
            dim_sq = size(mesh_occp_MarkerDen,2)
            if(present(mesh5D_norm_occp_dxdv)) then
               if(dim_sr/=size(mesh5D_norm_occp_dxdv,1)) then
                  write(stdout,*) "Error!! mesh5D_MarkerDen's dim_sr=",dim_sr," is different from mesh5D_dxdv's dim_s=",size(mesh5D_norm_occp_dxdv,1)
                  call flush(stdout)
                  stop
               endif
               if(dim_sq/=size(mesh5D_norm_occp_dxdv,2)) then
                  write(stdout,*) "Error!! mesh5D_MarkerDen's dim_sq=",dim_sq," is different from mesh5D_dxdv's dim_q=",size(mesh5D_norm_occp_dxdv,2)
                  call flush(stdout)
                  stop
               endif
            endif
         elseif(present(mesh5D_norm_occp_dxdv)) then
            dim_sr = size(mesh5D_norm_occp_dxdv,1)
            dim_sq = size(mesh5D_norm_occp_dxdv,2)
         endif


         ! if(mype==0) then
         !    write(stdout,'(2(A,I0))') "In Generate_Electrons_in_given_Phase_Space_PolarVelCoord_AdapTgtUr, dim_sr=",dim_sr,"; dim_sq=",dim_sq
         !    call flush(stdout)
         ! endif

         a_min = minval(a); a_max = maxval(a); dim_a = size(a)
         gq_min = minval(gq); gq_max = maxval(gq); dim_gq = size(gq)
         gf_min = minval(gf); gf_max = maxval(gf); dim_gf = size(gf)

         Del_sr = sr_max-sr_min;    dsr=Del_sr/(dim_sr-1)
         Del_sq = sq_max-sq_min;    dsq=Del_sq/(dim_sq-1)
         Del_a = a_max-a_min;       da=Del_a/(dim_a-1)
         Del_gq = gq_max-gq_min;    dgq=Del_gq/(dim_gq-1)
         Del_gf = gf_max-gf_min;    dgf=Del_gf/(dim_gf-1)


         allocate(sr_1D(dim_sr),sq_1D(dim_sq))

         do srid=1,dim_sr
            sr_1D(srid)=sr_min+(srid-1)*dsr
         enddo
         do sqid=1,dim_sq
            sq_1D(sqid)=sq_min+(sqid-1)*dsq
         enddo

         ini_MarkerDen_mype = ini_me/((a1-a0)*(2*pi)*(zetamax-zetamin)*0.5*pi*umax**2.0)
         ini_MarkerDen_npartd = ini_me_npartd/((a1-a0)*(2*pi)*(zetamax-zetamin)*0.5*pi*umax**2.0)



         inv_da = 1.0/da; inv_dgq=1.0/dgq; inv_dsr=1.0/dsr; inv_dsq=1.0/dsq

         ! BLK_SZa = PtlSrcCV%dim_a-1
         ! BLK_SZq = PtlSrcCV%dim_gq-1
         BLK_SZa = 7
         BLK_SZq = 20
         ! BLK_SZa = 5
         ! BLK_SZq = 5
         num_BLKs = dim_gf*ceiling(real(dim_a-1)/BLK_SZa)*ceiling(real(dim_gq-1)/BLK_SZq)
         allocate(BLKs_aid(num_BLKs),BLKs_gqid(num_BLKs),BLKs_gfid(num_BLKs),BLKs_m0(num_BLKs), BLKs_m1(num_BLKs))
         allocate(BLKs_Nsamp_mype(num_BLKs), BLKs_Nsamp_npartd(num_BLKs),BLKs_NnewP_mype(num_BLKs),BLKs_NnewP_npartd(num_BLKs))
         allocate(BLKs_d3x(num_BLKs),BLKs_d2u(num_BLKs))


         Bid=0
         do fid=1,dim_gf
            do aid=1,dim_a-1,BLK_SZa
               do qid=1,dim_gq-1,BLK_SZq
                  Bid = Bid + 1
                  BLKs_aid(Bid)=aid
                  BLKs_gqid(Bid)=qid
                  BLKs_gfid(Bid)=fid
               enddo
            enddo
         enddo

         if(Bid/=num_BLKs) then
            write(stdout,'((A,I0,A,I0))') "Error!! Something is Wrong!! After finding block idices, Bid=", Bid, "; num_BLKs=", num_BLKs
            stop
         endif

         m0 = me+1

         do Bid=1,num_BLKs

            BLKs_m0(Bid) = me+1

            aid = BLKs_aid(Bid)
            qid = BLKs_gqid(Bid)
            fid = BLKs_gfid(Bid)

            loc_max_aid = min(dim_a,aid+BLK_SZa)
            loc_Del_a = a(loc_max_aid)-a(aid)
            loc_max_qid = min(dim_gq,qid+BLK_SZq)
            loc_Del_q = gq(loc_max_qid)-gq(qid)

            loc_gf_max = min(gf(dim_gf),gf(fid)+0.5*dgf)
            loc_gf_min = max(gf(1),gf(fid)-0.5*dgf)
            loc_Del_f = loc_gf_max - loc_gf_min

            BLKs_d3x(Bid) = loc_Del_f*loc_Del_a*loc_Del_q

            loc_min_ur = minval(tgt_ur_bound(aid:loc_max_aid, qid:loc_max_qid, fid,1)) ! lower bound
            loc_max_ur = maxval(tgt_ur_bound(aid:loc_max_aid, qid:loc_max_qid, fid,2)) ! upper bound

            BLKs_d2u(Bid) = 0.5*pi*(loc_max_ur**2.0-loc_min_ur**2.0)

            ! loc_Nsamp = ceiling(ini_Ne_marker_density*loc_Del_a*loc_Del_q*loc_Del_f*0.5*pi*(loc_max_ur**2.0-loc_min_ur**2.0) +sprng())
            loc_Nsamp = ceiling(ini_MarkerDen_mype*BLKs_d3x(Bid)*BLKs_d2u(Bid) +sprng())
            BLKs_Nsamp_mype(Bid) = loc_Nsamp

               ! write(stdout,'((A,I0,A,3(I0,1x),A,I0,2(A,ES9.3)))') "mype=",mype,"; aid,qid,fid=", aid, qid,fid,"; Bid=",Bid, "; loc_min_ur=", loc_min_ur, "; loc_max_ur=",loc_max_ur
               ! call flush(stdout)


            ! if(mype==0) then
            !    write(stdout,'(4(A,I0),2(A,ES9.3))') "num_BLKs=",num_BLKs, "; Bid=", Bid, "; loc_Nsamp=",loc_Nsamp,"; me=", me,"; loc_max_ur=", loc_max_ur, "; loc_min_ur=", loc_min_ur
            !    call flush(stdout)
            ! endif

            BLKs_NnewP_mype(Bid) = 0
            do k=1,loc_Nsamp
               ptl_a  = a(aid) +loc_Del_a*sprng()
               ptl_gq = gq(qid)+loc_Del_q*sprng() ! NOTE!! This is gq value at either zeta=zetamin or zeta=zetamax

               ! ptl_aid  = aid+floor((ptl_a-a(aid))*inv_da+0.5) ! Find nearst grid
               ! ptl_gqid = qid+floor((ptl_gq-gq(qid))*inv_dgq+0.5) ! Find nearst grid

               ptl_aid  = ceiling((ptl_a-a_min)*inv_da+0.5) ! Find nearst grid
               ptl_gqid = ceiling((ptl_gq-gq_min)*inv_dgq+0.5) ! Find nearst grid

               ! ur lower(1) and upper(2) bound at ptl position
               ptl_ur_bound = tgt_ur_bound(ptl_aid,ptl_gqid,fid,:)

               ! Find sampling points uniformly distributed within an annulus in the poloar coordinates (vpara,vperp)
               ! r=sqrt(a1^2+rand*(a2^2-a1^2));
               ptl_ur = sqrt(loc_min_ur**2.0 + sprng()*(loc_max_ur**2.0-loc_min_ur**2.0))

               ! if(mype==0) then
               !    write(stdout,'((A,I0,A,3(I0,1x),3(A,I0),A,ES8.2,A,2(ES8.2,1x),2(A,ES9.2)))') "mype=",mype,"; aid,qid,fid=", aid, qid,fid,"; Bid=",Bid, "; loc_Nsamp=",loc_Nsamp,"; ptl_aid=",ptl_aid, "; ptl_a=", ptl_a,"; ptl_ur_bound=",ptl_ur_bound(:),"; loc_min_ur=", loc_min_ur, "; loc_max_ur=",loc_max_ur
               !    call flush(stdout)
               ! endif


               if( (ptl_ur-ptl_ur_bound(1))*(ptl_ur-ptl_ur_bound(2)) < 0 ) then
                  ! sr (speed) refers normalized ur regarding tgt_ur_bound:
                  ! sr(sr,sq,a,gq,gf) = (ur(sr,sq,a,gq,gf)-ur_LowerBound(a,gq,gf))/(ur_UpperBound(a,gq,gf)-ur_LowerBound(a,gq,gf))
                  ptl_sr = (ptl_ur-ptl_ur_bound(1))/(ptl_ur_bound(2)-ptl_ur_bound(1))
                  ptl_srid = ceiling((ptl_sr-sr_min)*inv_dsr+0.5)

                  ptl_sq = sq_min + sprng()*(sq_max-sq_min)
                  ptl_sqid = ceiling((ptl_sq-sq_min)*inv_dsq+0.5)

                  !accept this particle
                  BLKs_NnewP_mype(Bid) = BLKs_NnewP_mype(Bid) +1
                  me = me+1

                  ! This is only a part of zeon0(6,m); it will be recalulated by multiplying 1/ge_val later
                  if(present(mesh5D_norm_occp_dxdv)) then
                     zeon0(6,me) = 1.0 -mesh5D_norm_occp_dxdv(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid)
                  else
                     zeon0(6,me) = 1.0
                  endif

                  if(present(mesh_occp_MarkerDen)) then

                     sq0 = max(sq_min, sq_1D(ptl_sqid)-0.5*dsq)
                     sq1 = min(sq_max, sq_1D(ptl_sqid)+0.5*dsq)

                     sr0 = max(sr_min, sr_1D(ptl_srid)-0.5*dsr)
                     sr1 = min(sr_max, sr_1D(ptl_srid)+0.5*dsr)
                     ur0 = ptl_ur_bound(1) + sr0*(ptl_ur_bound(2)-ptl_ur_bound(1))
                     ur1 = ptl_ur_bound(1) + sr1*(ptl_ur_bound(2)-ptl_ur_bound(1))

                     ! ! mesh_cell_area =  Cell_dVol_x3*0.5*(ur1**2.0-ur0**2.0)*(sq1-sq0)
                     tgt_cell_Nm = ini_MarkerDen_mype*(min(a1,a(ptl_aid)+0.5*da)-max(a0,a(ptl_aid)-0.5*da))*dgq*(dgf)*0.5*(ur1**2.0-ur0**2.0)*(sq1-sq0)

                     ! if(mype==0) then
                     !    if(ur1<umax) then
                     !       write(stdout,'(2(A,ES10.4))') "tgt_cell_Nm=", tgt_cell_Nm, "; mesh_MarkerDen_mype=",mesh_occp_MarkerDen(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid)
                     !       call flush(stdout)
                     !    endif
                     ! endif

                     ! if(mype==0) then
                     !    if(ptl_aid>0.85 ) then
                     !       write(stdout,'(3(A,ES10.4))') "tgt_cell_Nm=", tgt_cell_Nm, "; ini_MarkerDen_npartd=",ini_MarkerDen_npartd,"; mesh_MarkerDen=",mesh_occp_MarkerDen(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid)
                     !       call flush(stdout)
                     !    endif
                     ! endif

                     ! if(tgt_cell_Nm<1e-3) then
                     if(.false.) then
                        ! if tgt_cell_Nm is too small, it needs some amplification for a more proper determination.
                        ! amp_factor =(min(dim_sr,ptl_srid+1)-max(1,ptl_srid-1)+1)*(min(dim_sq,ptl_sqid+1)-max(1,ptl_sqid-1)+1)*(min(dim_a,ptl_aid+3)-max(1,ptl_aid-3))*(min(dim_gq,ptl_gqid+6)- max(1,ptl_gqid-6))*dim_gf
                        ! amp_tgt_cell_Nm = amp_factor*tgt_cell_Nm
                        ! amp_occp_MarkerNum_mype = (1.0/npartdom)*sum(mesh_occp_MarkerNum_npartd(max(1,ptl_srid-1):min(dim_sr,ptl_srid+1), max(1,ptl_sqid-1):min(dim_sq,ptl_sqid+1), max(1,ptl_aid-3):min(dim_a,ptl_aid+3), max(1,ptl_gqid-6):min(dim_gq,ptl_gqid+6), :))
                        ! amp_factor =(min(dim_a,ptl_aid+3)-max(1,ptl_aid-3)+1)*(min(dim_gq,ptl_gqid+6)- max(1,ptl_gqid-6)+1)*dim_gf
                        ! amp_tgt_cell_Nm = c*tgt_cell_Nm
                        ! amp_occp_MarkerNum_mype = (1.0/npartdom)*sum(mesh_occp_MarkerDen(ptl_srid, ptl_sqid, max(1,ptl_aid-3):min(dim_a,ptl_aid+3), max(1,ptl_gqid-6):min(dim_gq,ptl_gqid+6), :))
                        ! amp_factor =(min(dim_a,ptl_aid+3)-max(1,ptl_aid-3)+1)*(min(dim_gq,ptl_gqid+6)- max(1,ptl_gqid-6)+1)
                        ! mean_occp_MarkerDen = (1.0/amp_factor)*sum(mesh_occp_MarkerDen(ptl_srid, ptl_sqid, max(1,ptl_aid-3):min(dim_a,ptl_aid+3), max(1,ptl_gqid-6):min(dim_gq,ptl_gqid+6), fid))

                        amp_factor =(min(dim_a,ptl_aid+5)-max(1,ptl_aid-5)+1)*(min(dim_gq,ptl_gqid+10)- max(1,ptl_gqid-10)+1)
                        mean_occp_MarkerDen = (1.0/amp_factor)*sum(mesh_occp_MarkerDen(ptl_srid, ptl_sqid, max(1,ptl_aid-5):min(dim_a,ptl_aid+5), max(1,ptl_gqid-10):min(dim_gq,ptl_gqid+10), fid))

                        if(present(mesh5D_norm_occp_dxdv)) then
                           ! mean_norm_occp_dxdv = (1.0/amp_factor)*sum(mesh5D_norm_occp_dxdv(ptl_srid, ptl_sqid, max(1,ptl_aid-3):min(dim_a,ptl_aid+3), max(1,ptl_gqid-6):min(dim_gq,ptl_gqid+6), fid))
                           mean_norm_occp_dxdv = (1.0/amp_factor)*sum(mesh5D_norm_occp_dxdv(ptl_srid, ptl_sqid, max(1,ptl_aid-5):min(dim_a,ptl_aid+5), max(1,ptl_gqid-10):min(dim_gq,ptl_gqid+10), fid))
                           zeon0(6,me) = 1.0 -mean_norm_occp_dxdv
                        endif

                        if(mype==0) then
                           if(ptl_a>0.85) then
                              write(stdout,'(4(A,ES9.2))') "tgt_cell_Nm=", tgt_cell_Nm,"; amp_tgt_cell_Nm=", amp_factor*tgt_cell_Nm, "; dxdv=",mesh5D_norm_occp_dxdv(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid),"; mean_dxdv=", mean_norm_occp_dxdv
                              call flush(stdout)
                           endif
                        endif


                        ! if(mype==0) then
                        !    if(ur1<umax) then
                        !       write(stdout,'(2(A,ES10.4))') "amp_factor*tgt_cell_Nm=",amp_factor*tgt_cell_Nm,"; ini_MarkerDen_npartd=", ini_MarkerDen_npartd, "; mean_occp_MarkerNum=",mean_occp_MarkerNum
                        !       call flush(stdout)
                        !    endif
                        ! endif

                        if(mean_occp_MarkerDen>0.95*ini_MarkerDen_npartd .and. mean_occp_MarkerDen<=ini_MarkerDen_npartd) then
                           mean_occp_MarkerDen=0.95*ini_MarkerDen_npartd ! To prevent too close occp_MarkerDen~ini_MarkerDen which could be troublesome
                        endif

                        if(sprng()*ini_MarkerDen_npartd < mean_occp_MarkerDen) then
                           BLKs_NnewP_mype(Bid) = BLKs_NnewP_mype(Bid) -1
                           me=me-1
                           cycle ! reject this particle
                        endif
                        ! Note that denominator is always smaller than 1 otherwise it should be rejected by the previous checking routine
                        zeon0(6,me) = zeon0(6,me)/(1.0-mean_occp_MarkerDen/ini_MarkerDen_npartd)



                     else
                        ! occp_MarkerDen = (1.0/npartdom)*mesh_occp_MarkerDen(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid)
                        occp_MarkerDen = mesh_occp_MarkerDen(ptl_srid, ptl_sqid, ptl_aid, ptl_gqid, fid)

                        if(occp_MarkerDen>0.95*ini_MarkerDen_npartd .and. occp_MarkerDen<=ini_MarkerDen_npartd) then
                           occp_MarkerDen=0.95*ini_MarkerDen_npartd ! To prevent too close occp_MarkerDen~ini_MarkerDen which could be troublesome
                        endif

                        if(sprng()*ini_MarkerDen_npartd <occp_MarkerDen) then
                           BLKs_NnewP_mype(Bid) = BLKs_NnewP_mype(Bid) -1
                           me=me-1
                           cycle ! reject this particle
                        endif
                        ! Note that denominator is always smaller than 1 otherwise it should be rejected by the previous checking routine
                        zeon0(6,me) = zeon0(6,me)/(1.0-occp_MarkerDen/ini_MarkerDen_npartd)
                     endif
                  endif


                  zeon(1,me) = ptl_a
                  zeon(7,me) = D01AHF(0.0_wp,zeon(1,me),error_8,minpts,ac,dsda,-1,ifail)

                  ! Interpolate safety factor
                  tmp_aid = max(1,min(int((zeon(1,me)-PSCV%a_min)*PSCV%inv_da) +1,PSCV%dim_a-1)) ! +1 is due to Fortran's default index starting from 1
                  ws1 = (zeon(1,me)-PSCV%a_1D(tmp_aid))*PSCV%inv_da
                  ws0 = (1.0-ws1)
                  P_safety_q = ws0*PSCV%safety_q_1D(tmp_aid) + ws1*PSCV%safety_q_1D(tmp_aid+1)

                  zeon(3,me) = modulo(loc_gf_min + loc_Del_f*sprng(), 2.0*pi)
                  ! zeon(2,me) = modulo(ptl_gq + (gf(fid)-zeon(3,me))/P_safety_q, 2.0*pi)
                  zeon(2,me) = modulo(ptl_gq + (zeon(3,me)-gf(fid))/P_safety_q, 2.0*pi)

                  upara = (ptl_ur*cos(ptl_sq))
                  uperp = (ptl_ur*sin(ptl_sq))

                  ! Note that zeon(4) and zeon(6) are temporary used to carray upara and uperp values
                  ! The real zeon values will be calculated later
                  zeon(4,me) = upara
                  zeon(6,me) = uperp

                  elabels(me) = label_val
               endif
            enddo ! End k=1,loc_Nsamp

            BLKs_m1(Bid) = me

            ! if(mype==0) then
            !    write(stdout,'(4(A,I0),A,F8.3)') "num_BLKs=",num_BLKs, "; Bid=", Bid, "; this_BLK_NnewP_mype=",BLKs_NnewP_mype(Bid),"; this_BLKs_Nsamp_npartd=",BLKs_Nsamp_mype(Bid),"; Acceptance Rate=",real(BLKs_NnewP_mype(Bid))/(BLKs_Nsamp_mype(Bid))
            !    call flush(stdout)
            ! endif

         enddo ! End BLK loop
         !$acc update device(zeon(:,m0:me),elabels(m0:me))
         call MPI_ALLREDUCE(BLKs_NnewP_mype, BLKs_NnewP_npartd, num_BLKs, MPI_INTEGER, MPI_SUM, partd_comm,ierror)
         call MPI_ALLREDUCE(BLKs_Nsamp_mype, BLKs_Nsamp_npartd, num_BLKs, MPI_INTEGER, MPI_SUM, partd_comm,ierror)

         if(mype==0) then
            write(stdout,'(3(A,I0),A,F8.3)') "num_BLKs=",num_BLKs, "; tot_BLKs_NnewP_npartd=",sum(BLKs_NnewP_npartd),"; tot_BLKs_Nsamp_npartd=",sum(BLKs_Nsamp_npartd),"; Acceptance Rate=",real(sum(BLKs_NnewP_npartd))/sum(BLKs_Nsamp_npartd)
            write(stdout,'(3(A,I0))') "tot_BLKs_NnewP_mype=",sum(BLKs_NnewP_mype),"; me=",me,"; memax=",memax
            call flush(stdout)
         endif

         call n_T_u_group_e(m0,me)   ! get Te, note that back_f0_e cannot be calculated here because it requires zeon(4) and zeon(6)
         !$acc update self(backTe(m0:me))
         call esigetbf_partial(zeon,m0,me,zeondim)

         do Bid=1,num_BLKs
            do m=BLKs_m0(Bid),BLKs_m1(Bid)
               upara = zeon(4,m) ! upara value is carried by zeon(4) from the previous do loop
               uperp = zeon(6,m) ! uperp value is carried by zeon(6) from the previous do loop

               vpara = upara*sqrt(backTe(m))
               vperp = uperp*sqrt(backTe(m))

               zeon(4,m) = msm0*vpara/(csc0*BB(m))
               zeon(6,m) = vperp*sqrt(msm0/(2.0*BB(m)))

               ! Note that ge_val (marker distribution) should include the Jacobian from the polar coordinate (here, vperp)
               ! ge_val(m) = tot_N_newP/(Del_a*Del_gq*Del_gf*(0.5*pi*(ur_max**2.0-ur_min**2.0))*backTe(m))/(Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)
               ! ge_val = BLKs_NnewP_npartd(Bid)/(BLKs_d3x(Bid)*BLKs_d2u(Bid)*backTe(m)*Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)
               ge_val = BLKs_Nsamp_npartd(Bid)/(BLKs_d3x(Bid)*BLKs_d2u(Bid)*backTe(m)*Jacobian(zeon(1,m),zeon(2,m))*2.0*pi*vperp)
               zeon0(6,m) = zeon0(6,m)/ge_val ! Phase space volume
            enddo
         enddo
         !$acc update device(zeon(:,m0:me))
         !$acc update device(zeon0(6,m0:me))
         call n_T_u_group_e(m0,me,flag_backf0=.true.)   ! This is to get back_f0_e inside. (which requires zeon(4) and zeon(6))


         ! if (mype==0) then
         !    write(filename,'("test_",i5.5,".dat")') istep
         !    open(unit=20211210, file=filename, position='append',action='write')
         !    do m=m0,me
         !       upara = zeon(4,m)*csc0*BB(m)*sqrt(m0ms/(msm0*backTe(m)))! upara
         !       uperp = zeon(6,m)*sqrt(2.0*BB(m))/sqrt(msm0*backTe(m))  ! uperp
         !       write(20211210,'(99(ES16.8,1x))') zeon(:,m), upara, uperp, zeon0(6,m)
         !    enddo
         !    call flush(20211210)
         !    close(20211210)
         ! endif




         ! ! Calculate Adiab term
         ! allocate(AdiabTerm(m0:me))
         ! call electron_locator_partial(m0,me) ! This is required to calculate wpeon, wteon0 and etc
         ! !$acc parallel loop present(AdiabTerm, wzeon,jteon0,wpeon,wteon0,jteon1,wteon1)
         ! do m=m0,me
         !    wz1=wzeon(m)
         !    wz0=1.0_wp-wz1

         !    ij=jteon0(m)
         !    wp0=1.0_wp-wpeon(m)
         !    wt00=1.0_wp-wteon0(m)
         !    intp_phi = wp0*wt00*(wz0*phi(0,ij)+wz1*phi(1,ij))
         !    intp_FLA_Te0 = wp0*wt00*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

         !    ij=ij+1
         !    wt10=1.0_wp-wt00
         !    intp_phi=intp_phi+wp0*wt10*(wz0*phi(0,ij)+wz1*phi(1,ij))
         !    intp_FLA_Te0=intp_FLA_Te0+wp0*wt10*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

         !    ij=jteon1(m)
         !    wp1=1.0_wp-wp0
         !    wt01=1.0_wp-wteon1(m)
         !    intp_phi=intp_phi+wp1*wt01*(wz0*phi(0,ij)+wz1*phi(1,ij))
         !    intp_FLA_Te0=intp_FLA_Te0+wp1*wt01*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

         !    ij=ij+1
         !    wt11=1.0_wp-wt01
         !    intp_phi=intp_phi+wp1*wt11*(wz0*phi(0,ij)+wz1*phi(1,ij)) ! --> phi at electrons
         !    intp_FLA_Te0=intp_FLA_Te0+wp1*wt11*(wz0*FLA_Te0(0,ij)+wz1*FLA_Te0(1,ij))

         !    r=(zeon(7,m)-s0)/deltar       ! in radial coordinator s
         !    ip=max(0,min(mpsi-1,int(r)))      ! ip [0,mpsi-1]
         !    fip=r-real(ip,wp)
         !    intp_phi00=fip*phi00(ip+1)+(1.0_wp-fip)*phi00(ip) ! <phi> at electrons
         !    if ( Flag_Def_dne_adiab==0 ) then
         !       AdiabTerm(m)=abs(csc0)*(intp_phi-real(1.0-ZF_SEPARATION,wp)*intp_phi00)/(backTe(m)*msm0) ! Adiabatic term (e*dphi_NZ/Te)
         !    elseif(Flag_Def_dne_adiab==1) then
         !       ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) for the New Poisson3D solver
         !       AdiabTerm(m)=abs(csc0)*(intp_phi+real(ZF_SEPARATION,wp)*intp_phi00)/(backTe(m)*msm0) ! Adiabatic term (e*dphi_tot/Te)
         !    elseif(Flag_Def_dne_adiab==2) then
         !       ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) with FLA_Te0
         !       AdiabTerm(m)=abs(csc0)*(intp_phi+real(ZF_SEPARATION,wp)*intp_phi00)/(intp_FLA_Te0) ! Adiabatic term (e*dphi_tot/Te)
         !    elseif(Flag_Def_dne_adiab==3) then
         !       ! MGYOO (02/19/2020): New Adiabatic term using total phi (zonal+nonzonal) with FLA_Te0
         !       AdiabTerm(m)=-1.0+exp(abs(csc0)*(intp_phi+real(ZF_SEPARATION,wp)*intp_phi00)/(backTe(m)*msm0)) ! Adiabatic term (e*dphi_tot/Te)
         !    else
         !       AdiabTerm(m)=-1e100
         !       write(stdout,*) "Error!! Flag_Def_dne_adiab=",Flag_Def_dne_adiab, " is not defined in pushe_new.F90!!"
         !       STOP
         !    endif
         ! enddo




         ! Set initial weights of ions
         ! zeon(5,m) is "df"
         ! zeon(8,m) is "f_tot"
         ! zeon0(6,m) is "phase volume"

         !$acc parallel loop present(zeon,zeon0,back_f0_e,we_init,we_vpara,we_vperp,we_accel,we_fdamp,UnNorm_theta_e,UnNorm_zeta_e)
         do m=m0,me
            ! zeon(5,m) is "df" or "dh"
            ! zeon(8,m) is "f_tot"
            ! zeon(9,m) is "f_AD"
            ! zeon0(9,m) is "f_tot"
            ! zeon0(6,m) is "phase volume"

            integral_f0 = zeon0(6,m)*back_f0_e(m) !f_0

            target_w_tot = integral_f0



            ! zeon(9,m) is f_AD(@ irk=1) which is the adiabatic response (=f0*e*phi/Te0) or (=f0*[exp(e*phi/Te0)-1])
            ! zeon(9,m) = integral_f0*AdiabTerm(m) ! f_AD


            call Get_Intp_phi_on_SG(zeon(7,m), zeon(2,m), zeon(3,m), intp_phi_V)
            Te_inverse=1.0_wp/backTe(m)*m0ms  ! normalized-Te=Te/E_0, the same unit as normalized Ti
            AD = (Fac_Volt2GTS*intp_phi_V)*Te_inverse*abs(csc0)
            zeon(9,m) = integral_f0*AD ! f_AD



            zeon(8,m) = target_w_tot -int(ifelectron/2)*zeon(9,m)! full-f (ifelectron<2) or full-h (ifelectron==2)
            zeon(5,m)  = zeon(8,m)-integral_f0 ! df or dh
            zeon0(9,m) = target_w_tot ! f_tot

            we_init(m)  = zeon(5,m)
            we_vpara(m) = 0
            we_vperp(m) = 0
            we_accel(m) = 0
            we_adiab(m) = 0
            we_fdamp(m) = 0

            !$acc loop seq
            do i=1,5
               zeon0(i,m)     = zeon(i,m)
            enddo
            zeon0(7,m)       = zeon(7,m)
            zeon0(8,m)       = zeon(8,m)
            UnNorm_theta_e(m)= zeon(2,m)
            UnNorm_zeta_e(m) = zeon(3,m)
         enddo
         acc_newer_D_zeon=.true.


         ! deallocate(AdiabTerm)

         deallocate(BLKs_aid,BLKs_gqid,BLKs_gfid,BLKs_m0,BLKs_m1)
         deallocate(BLKs_Nsamp_mype, BLKs_Nsamp_npartd,BLKs_NnewP_mype,BLKs_NnewP_npartd)
         deallocate(BLKs_d3x,BLKs_d2u)
         deallocate(sr_1D,sq_1D)


#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Generate_Electrons_in_given_Phase_Space_PolarVelCoord_AdapTgtUr takes", time1-time0, " secs"
         call flush(stdout)
#endif

      end subroutine Generate_Electrons_in_given_Phase_Space_PolarVelCoord_AdapTgtUr

      subroutine Clear_newIons_beyond_Velocity_Boundary(m_start, m_end)
         use global_parameters
         use particle_array
         use esi_interface
         use function_gg_module
         use mpi_defs
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0, m1
         integer dim_zion, dim_zion0, dim_RK_cum_zion
         integer, allocatable :: clear_mid_arr(:), fill_mid_arr(:)

         real(wp), allocatable :: zion_4_arr(:),zion_6_arr(:) ! temporary array 4 (rho_para) and 6 (sqrt(mu))
         real(wp) upara, uperp
         integer N_clear, i,j,m,mtop,lasth,fill_num
         real(wp) time0,time1


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

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

         call n_T_u_group(m0,m1)   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(backu(m0:m1),backT(m0:m1))
         call esigetbf_partial(zion,m0,m1,ziondim)
         backu(m0:m1)=backu(m0:m1)*F(m0:m1)/BB(m0:m1) ! =u_//=I*omega_t/B

         ! allocate(zion_4_arr(m0:m1),zion_6_arr(m0:m1))
         allocate(clear_mid_arr(m1-m0+1))

         ! Update ilabels and zion(4,:), zion(6,:) from device to host
         !$acc update self(ilabels(m0:m1),zion(:,m0:m1))


         ! !$acc parallel loop
         ! do m=m0,m1
         !    zion_4_arr(m)=zion(4,m)
         !    zion_6_arr(m)=zion(6,m)
         ! enddo


         ! Find indices of target ions to clear
         N_clear=0
         do m=m0,m1
            if(ilabels(m)<0) then
               ! upara = (zion_4_arr(m)*BB(m) -backu(m))/sqrt(backT(m))
               ! uperp = zion_6_arr(m)*sqrt(2.0*BB(m)/backT(m))
               upara = (zion(4,m)*BB(m) -backu(m))/sqrt(backT(m))
               uperp = zion(6,m)*sqrt(2.0*BB(m)/backT(m))
               if(Flag_Load_and_WeightDef==1 .or. Flag_Load_and_WeightDef==2) then
                  if(sqrt(upara**2.0+uperp**2.0)>umax) then
                     N_clear=N_clear+1
                     clear_mid_arr(N_clear) = m
                  endif
               endif
            endif
         enddo


         mtop = mi
         mi = mi-N_clear
         lasth = N_clear

         allocate(fill_mid_arr(N_clear))
         fill_num=0
         do i=1,N_clear
            m = clear_mid_arr(i)
            if(m>mi) exit
            do while(mtop == clear_mid_arr(lasth))
               mtop=mtop-1
               lasth=lasth-1
            enddo
            fill_mid_arr(i) = mtop
            mtop = mtop-1
            fill_num = i
         enddo


         !$acc parallel loop
         do i=1, fill_num
            m = clear_mid_arr(i)
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
         acc_newer_D_zion=.true.

         call ion_locator_at_larmor_points_partial(m_start=m0, m_end=mi) ! This may not be necessary, but just in case something goes wrong in the future
         call n_T_u_group(m_start=m0, m_end=mi)

         ! deallocate(zion_4_arr,zion_6_arr)
         deallocate(clear_mid_arr)
         deallocate(fill_mid_arr)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Clear_newIons_beyond_Velocity_Boundary takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Clear_newIons_beyond_Velocity_Boundary

      subroutine Clear_some_Ions_beyond_Boundaries(m_start, m_end)
         use global_parameters
         use particle_array
         use esi_interface
         use function_gg_module
         use mpi_defs
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0, m1
         integer dim_zion, dim_zion0, dim_RK_cum_zion
         integer, allocatable :: clear_mid_arr(:), fill_mid_arr(:)

         real(wp) P_a
         real(wp) upara, uperp
         integer N_clear, i,j,m,mtop,lasth,fill_num
         real(wp) time0,time1


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

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

         call n_T_u_group(m0,m1)   ! get Te, note that back_f0_i cannot be calculated here because it requires zion(4) and zion(6)
         !$acc update self(backu(m0:m1),backT(m0:m1))
         call esigetbf_partial(zion,m0,m1,ziondim)
         backu(m0:m1)=backu(m0:m1)*F(m0:m1)/BB(m0:m1) ! =u_//=I*omega_t/B

         allocate(clear_mid_arr(m1-m0+1))

         ! Update ilabels and zion(4,:), zion(6,:) from device to host
         !$acc update self(ilabels(m0:m1),zion(:,m0:m1))

         ! Find indices of target ions to clear
         N_clear=0
         do m=m0,m1
            P_a = zion(1,m)
            if( (P_a-a0)*(P_a-a1)>0) then
               N_clear=N_clear+1
               clear_mid_arr(N_clear) = m
            elseif(ilabels(m)<0) then
               upara = (zion(4,m)*BB(m) -backu(m))/sqrt(backT(m))
               uperp = zion(6,m)*sqrt(2.0*BB(m)/backT(m))
               if(Flag_Load_and_WeightDef==1 .or. Flag_Load_and_WeightDef==2) then
                  if(sqrt(upara**2.0+uperp**2.0)>umax) then
                     if(ilabels(m)==-2.0) then
                        N_clear=N_clear+1
                        clear_mid_arr(N_clear) = m
                     elseif(ilabels(m)==-1.0 .and. Ptl_BC0_val==0.0) then
                        ! when the particle comes from PSCV_BDY0 (left boundary)
                        ! and when f0(PSCV_BDY0)=0
                        N_clear=N_clear+1
                        clear_mid_arr(N_clear) = m
                     elseif(ilabels(m)==-3.0 .and. Ptl_BC1_val==0.0) then
                        ! when the particle comes from PSCV_BDY1 (right boundary)
                        ! and when f0(PSCV_BDY1)=0
                        N_clear=N_clear+1
                        clear_mid_arr(N_clear) = m
                     endif
                  endif
               endif
            endif
         enddo


         mtop = mi
         mi = mi-N_clear
         lasth = N_clear

         allocate(fill_mid_arr(N_clear))
         fill_num=0
         do i=1,N_clear
            m = clear_mid_arr(i)
            if(m>mi) exit
            do while(mtop == clear_mid_arr(lasth))
               mtop=mtop-1
               lasth=lasth-1
            enddo
            fill_mid_arr(i) = mtop
            mtop = mtop-1
            fill_num = i
         enddo


         !$acc parallel loop
         do i=1, fill_num
            m = clear_mid_arr(i)
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
         acc_newer_D_zion=.true.

         call ion_locator_at_larmor_points_partial(m_start=m0, m_end=mi) ! This may not be necessary, but just in case something goes wrong in the future
         call n_T_u_group(m_start=m0, m_end=mi)

         deallocate(clear_mid_arr)
         deallocate(fill_mid_arr)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Clear_newIons_beyond_Velocity_Boundary takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Clear_some_Ions_beyond_Boundaries

      subroutine Clear_newElectrons_beyond_Velocity_Boundary(m_start, m_end)
         use global_parameters
         use particle_array
         use esi_interface
         use function_gg_module
         use mpi_defs
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0, m1
         integer dim_zeon, dim_zeon0, dim_RK_cum_zeon
         integer, allocatable :: clear_mid_arr(:), fill_mid_arr(:)

         real(wp), allocatable :: zeon_4_arr(:),zeon_6_arr(:) ! temporary array 4 (rho_para) and 6 (sqrt(mu))
         real(wp) upara, uperp
         integer N_clear, i,j,m,mtop,lasth,fill_num
         real(wp) time0,time1


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

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

         call n_T_u_group_e(m0,m1)   ! get Te, note that back_f0_i cannot be calculated here because it requires zeon(4) and zeon(6)
         !$acc update self(backTe(m0:m1))
         call esigetbf_partial(zeon,m0,m1,zeondim) ! BB update to host

         allocate(zeon_4_arr(m0:m1),zeon_6_arr(m0:m1))
         allocate(clear_mid_arr(m1-m0+1))

         ! Update elabels and zeon(4,:), zeon(6,:) from device to host
         !$acc update self(elabels(m0:m1))
         !$acc parallel loop
         do m=m0,m1
            zeon_4_arr(m)=zeon(4,m)
            zeon_6_arr(m)=zeon(6,m)
         enddo


         ! Find indices of target ions to clear
         N_clear=0
         do m=m0,m1
            if(elabels(m)<0) then
               ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
               ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))
               upara = zeon_4_arr(m)*csc0*BB(m)*sqrt(m0ms/(msm0*backTe(m)))! upara
               uperp = zeon_6_arr(m)*sqrt(2.0*BB(m))/sqrt(msm0*backTe(m))  ! uperp

               if(Flag_Load_and_WeightDef==1 .or. Flag_Load_and_WeightDef==2) then
                  if(sqrt(upara**2.0+uperp**2.0)>umax) then
                     N_clear=N_clear+1
                     clear_mid_arr(N_clear) = m
                  endif
               endif
            endif
         enddo

         mtop = me
         me = me-N_clear
         lasth = N_clear

         allocate(fill_mid_arr(N_clear))
         fill_num=0
         do i=1,N_clear
            m = clear_mid_arr(i)
            if(m>me) exit
            do while(mtop == clear_mid_arr(lasth))
               mtop=mtop-1
               lasth=lasth-1
            enddo
            fill_mid_arr(i) = mtop
            mtop = mtop-1
            fill_num = i
         enddo


         !$acc parallel loop
         do i=1, fill_num
            m = clear_mid_arr(i)
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
         acc_newer_D_zeon=.true.

         call electron_locator_partial(m_start=m0, m_end=me) ! This may not be necessary, but just in case something goes wrong in the future
         call n_T_u_group_e(m_start=m0, m_end=me)

         deallocate(zeon_4_arr,zeon_6_arr)
         deallocate(clear_mid_arr)
         deallocate(fill_mid_arr)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Clear_newElectrons_beyond_Velocity_Boundary takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Clear_newElectrons_beyond_Velocity_Boundary

      subroutine Clear_some_Electrons_beyond_Boundaries(m_start, m_end)
         use global_parameters
         use particle_array
         use esi_interface
         use function_gg_module
         use mpi_defs
         implicit none

         integer, optional, intent(in) :: m_start, m_end
         integer m0, m1
         integer dim_zeon, dim_zeon0, dim_RK_cum_zeon
         integer, allocatable :: clear_mid_arr(:), fill_mid_arr(:)

         real(wp), allocatable :: zeon_146_arr(:,:) ! temporary array 1(a), 4 (rho_para) and 6 (sqrt(mu))
         real(wp) P_a, upara, uperp
         integer N_clear, i,j,m,mtop,lasth,fill_num
         real(wp) time0,time1


#ifdef _MGY_TIMER_
         time0=MPI_WTIME()
#endif

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

         call n_T_u_group_e(m0,m1)   ! get Te, note that back_f0_i cannot be calculated here because it requires zeon(4) and zeon(6)
         !$acc update self(backTe(m0:m1))
         call esigetbf_partial(zeon,m0,m1,zeondim) ! BB update to host

         allocate(zeon_146_arr(3,m0:m1))
         allocate(clear_mid_arr(m1-m0+1))

         ! Update elabels and zeon(4,:), zeon(6,:) from device to host
         !$acc update self(elabels(m0:m1))
         !$acc parallel loop
         do m=m0,m1
            zeon_146_arr(1,m)=zeon(1,m)
            zeon_146_arr(2,m)=zeon(4,m)
            zeon_146_arr(3,m)=zeon(6,m)
         enddo


         ! Find indices of target electrons to clear
         N_clear=0
         do m=m0,m1
            P_a = zeon_146_arr(1,m)
            if((P_a-a0)*(P_a-a1)>0) then
               N_clear=N_clear+1
               clear_mid_arr(N_clear) = m
            elseif(elabels(m)<0) then
               ! norm_vpara = vpara/sqrt(T/m) = rho_para*(eB/m)/sqrt(T/m) = rho_para*(e*B)/sqrt(m*T)
               ! norm_vperp = sqrt(mu)*sqrt(2B/m)/sqrt(T/m) = sqrt(mu)*sqrt(2B)/(sqrt(T))
               upara = zeon_146_arr(2,m)*csc0*BB(m)*sqrt(m0ms/(msm0*backTe(m)))! upara
               uperp = zeon_146_arr(3,m)*sqrt(2.0*BB(m))/sqrt(msm0*backTe(m))  ! uperp
               if(Flag_Load_and_WeightDef==1 .or. Flag_Load_and_WeightDef==2) then
                  if(sqrt(upara**2.0+uperp**2.0)>umax) then
                     if(elabels(m)==-2.0) then
                        N_clear=N_clear+1
                        clear_mid_arr(N_clear) = m
                     elseif(elabels(m)==-1.0 .and. Ptl_BC0_val==0.0) then
                        ! when the particle comes from PSCV_BDY0 (left boundary)
                        ! and when f0(PSCV_BDY0)=0
                        N_clear=N_clear+1
                        clear_mid_arr(N_clear) = m
                     elseif(elabels(m)==-3.0 .and. Ptl_BC1_val==0.0) then
                        ! when the particle comes from PSCV_BDY1 (right boundary)
                        ! and when f0(PSCV_BDY1)=0
                        N_clear=N_clear+1
                        clear_mid_arr(N_clear) = m
                     endif
                  endif
               endif
            endif
         enddo

         mtop = me
         me = me-N_clear
         lasth = N_clear

         allocate(fill_mid_arr(N_clear))
         fill_num=0
         do i=1,N_clear
            m = clear_mid_arr(i)
            if(m>me) exit
            do while(mtop == clear_mid_arr(lasth))
               mtop=mtop-1
               lasth=lasth-1
            enddo
            fill_mid_arr(i) = mtop
            mtop = mtop-1
            fill_num = i
         enddo


         !$acc parallel loop
         do i=1, fill_num
            m = clear_mid_arr(i)
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
         acc_newer_D_zeon=.true.

         call electron_locator_partial(m_start=m0, m_end=me) ! This may not be necessary, but just in case something goes wrong in the future
         call n_T_u_group_e(m_start=m0, m_end=me)

         deallocate(zeon_146_arr)
         deallocate(clear_mid_arr)
         deallocate(fill_mid_arr)

#ifdef _MGY_TIMER_
         time1=MPI_WTIME()
         if(mype==0) write(stdout,'(a,I5,a,ES11.3,a)') "mype=",mype, ", Clear_some_Electrons_beyond_Boundaries takes", time1-time0, " secs"
         call flush(stdout)
#endif
      end subroutine Clear_some_Electrons_beyond_Boundaries


end module Ptl_Src_Loss
