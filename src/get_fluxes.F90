subroutine get_fluxes(wpi,dgf,dgq)
   use precision
   use mpi_defs
   use global_parameters
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface

   implicit none

   real(wp)  wpi(3,mi),dgf(mi),dgq(mi)
   integer   m, i, j, k, kc, j1, j2, m1, m2, l1, l2, mc, k1, k2,ierr
   real(wp)  adot,thedot,cosx,sinx,aj,da_table,yy,part
   real(wp)  my_ion, total_ion

   real(wp), dimension(nn,nbin_a) :: ONE_me,TOTAL_pe
   !ZHIXIN: ONE_me_g extended from 7 to 9 dimension for particle, energy flux: fort.86,87
   !             dimension 10 added for parallel Reynolds stress 2016
   real(wp), dimension(10,0:mpsi) :: ONE_me_g,TOTAL_pe_g

   ! add for diagnosis of distribution functions
   real(wp), dimension(7,nprobe,-nv:nv,0:nv) :: df, df_tot
   ! df(1,:,:,:) --> particle distribution function
   ! df(2,:,:,:) --> distribution of toroidal momentum
   ! df(3,:,:,:) --> distribution of toroidal momentum flux
   ! df(4,:,:,:) --> distribution of energy
   ! df(5,:,:,:) --> distribution of energy flux
   ! df(6,:,:,:) --> distribution of particle flux
   ! df(7,:,:,:) --> distribution of bootstrap current
   real(wp), dimension(nprobe,-nv:nv,0:nv) :: df_heat   ! added 2009
   !
   real(wp) r_area,vmax,the2,the1,part_l1,the_area,vper,ss1,ss2,ss3,ss4,am,qq1,qq2,qq3,qq4,vpara0

   real(wp) factor,time,flux_p,flux_E,flux_q,flux_vz,flux_the,flux_PI,jb, &
      angular_m,diff_c,delta_sumv2,wmean_0,w2mean_0,w2mean_1, &
      vt,av2,Tg,dTi,avz,sum_n,sum_sn0,da_size,bootstrap,  &
      v2_lab, veloc, part_k1, rdum,  &
      En_tur_sum,En_00_sum,En_tot_sum,enk_sum,da_size_g,adum,TM_sum, &
      flux_PIp,flux_vp

   real(wp) flux_p_strm, flux_p_ExB, flux_E_strm, flux_E_ExB ! Radial Flux of particle and energy by streaming or ExB, respectively
   real(wp), external :: SMOOTHSTEP

   real(wp) r,vpara,energy

   real(wp) dpsi_inv,e_psi2,e_psi1,dpsi2_inv,dpsi1,dpsi2,dpsi1_inv

   factor=3000.0_wp*trans/pmass0

   da_size=(a1-a0)/real(nbin_a,wp)
   da_size_g=(s1-s0)/real(mpsi,wp)
   ONE_me=0.0_wp

   vmax=delta_v*real(nv)    ! maximum velocity in v_th unit
   df=0.0_wp
   df_tot=0.0_wp
   df_heat=0.0_wp
   my_ion=0.0_wp

   !!     call esigetb(zion0(1,:),zion0(2,:),mi)   ! get B
   !     call esigetb(zion0,mi)   ! get B

   !$omp parallel private(m,r,yy,kc,j1,j2,part,vpara,energy,vpara0) reduction(+:ONE_me,my_ion)
   !$omp do
   do m=1,mi

      if (Flag_w_Threshold==1 .and. (zion0(5,m)<w_Thold_min .or. zion0(5,m)>w_Thold_max)) cycle

      r=zion0(1,m)
      yy=(r-a0)/da_size
      kc=int(yy+0.5_wp)            !{0,nbin_a}
      kc=min(nbin_a,max(0,kc))      !{0,nbin_a}

      j1=max(1,kc)
      j2=min(nbin_a,kc+1)
      part=abs(real(j2,wp)-0.5_wp-yy)

      vpara=zion0(4,m)*BB(m)
      energy=wpi(2,m)                ! v^2/2

      ONE_me(1,j1)=ONE_me(1,j1)+part
      ONE_me(1,j2)=ONE_me(1,j2)+(1.0_wp-part)

      ONE_me(2,j1)=ONE_me(2,j1)+part*zion0(5,m)
      ONE_me(2,j2)=ONE_me(2,j2)+(1.0_wp-part)*zion0(5,m)

      my_ion=my_ion+zion0(5,m)

      ONE_me(3,j1)=ONE_me(3,j1)+part*energy
      ONE_me(3,j2)=ONE_me(3,j2)+(1.0_wp-part)*energy

      ONE_me(4,j1)=ONE_me(4,j1)+part*zion0(5,m)**2
      ONE_me(4,j2)=ONE_me(4,j2)+(1.0_wp-part)*zion0(5,m)**2

      ONE_me(5,j1)=ONE_me(5,j1)+part*zion0(5,m)*energy         ! Temp contained in delta-f
      ONE_me(5,j2)=ONE_me(5,j2)+(1.0_wp-part)*zion0(5,m)*energy

      ONE_me(6,j1)=ONE_me(6,j1)+part*zion0(5,m)*wpi(1,m)       ! particle flux; wpi(1,m)=dr/dt|_EXB, not da/dt|_EXB
      ONE_me(6,j2)=ONE_me(6,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)

      ONE_me(7,j1)=ONE_me(7,j1)+part*zion0(5,m)*wpi(1,m)*energy ! energy flux
      ONE_me(7,j2)=ONE_me(7,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*energy

      ONE_me(8,j1)=ONE_me(8,j1)+part*zion0(5,m)*vpara ! v_// in delta-f
      ONE_me(8,j2)=ONE_me(8,j2)+(1.0_wp-part)*zion0(5,m)*vpara

      !         ONE_me(9,j1)=ONE_me(9,j1)+part*zion0(5,m)*zion0(4,m)*F(m)  ! angular momentum
      !         ONE_me(9,j2)=ONE_me(9,j2)+(1.0_wp-part)*zion0(5,m)*zion0(4,m)*F(m)
      !
      ! replaced by:
      ONE_me(9,j1)=ONE_me(9,j1)+         part*zion0(5,m)*Pa(m)*Pa(m)*dgf(m)  ! angular momentum
      ONE_me(9,j2)=ONE_me(9,j2)+(1.0_wp-part)*zion0(5,m)*Pa(m)*Pa(m)*dgf(m)

      !        ONE_me(10,j1)=ONE_me(10,j1)+part*zion0(5,m)*wpi(1,m)*zion0(4,m)*F(m)        ! angular momentum flux
      !        ONE_me(10,j2)=ONE_me(10,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*zion0(4,m)*F(m)
      !
      ! replaced by:
      ONE_me(10,j1)=ONE_me(10,j1)+part*zion0(5,m)*wpi(1,m)*Pa(m)*Pa(m)*dgf(m)        ! angular momentum flux
      ONE_me(10,j2)=ONE_me(10,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*Pa(m)*Pa(m)*dgf(m)

      ONE_me(11,j1)=ONE_me(11,j1)+part*zion0(5,m)*vpara*BB(m) ! ion bootstrap current = v_//*B
      ONE_me(11,j2)=ONE_me(11,j2)+(1.0_wp-part)*zion0(5,m)*vpara*BB(m)

      !        ONE_me(12,j1)=ONE_me(12,j1)+part*zion0(5,m)*wpi(1,m)*TT(m)*dgq(m)              ! poloidal angular momentum flux; added Aug. 2012 with TT(m)=1/|grad(theta)|
      !        ONE_me(12,j2)=ONE_me(12,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*TT(m)*dgq(m)

      ONE_me(12,j1)=ONE_me(12,j1)+part*zion0(5,m)*wpi(1,m)*TT(m)              ! poloidal angular momentum flux; modified Oct. 2012 with TT(m)= poloidal velocity
      ONE_me(12,j2)=ONE_me(12,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*TT(m)

      ONE_me(13,j1)=ONE_me(13,j1)+part*zion0(5,m)*TT(m)              ! poloidal flow velocity; modified Oct. 2012 with TT(m)= poloidal velocity
      ONE_me(13,j2)=ONE_me(13,j2)+(1.0_wp-part)*zion0(5,m)*TT(m)

      ! MGYOO(08/07/2019)
      ONE_me(14,j1)=ONE_me(14,j1)+part*zion0(5,m)*vpara_dvi_a(m)       ! particle radial flux due to streaming
      ONE_me(14,j2)=ONE_me(14,j2)+(1.0_wp-part)*zion0(5,m)*vpara_dvi_a(m)

      ONE_me(15,j1)=ONE_me(15,j1)+part*zion0(5,m)*vpara_dvi_a(m)*energy ! energy radial flux due to streaming
      ONE_me(15,j2)=ONE_me(15,j2)+(1.0_wp-part)*zion0(5,m)*vpara_dvi_a(m)*energy

      ONE_me(16,j1)=ONE_me(16,j1)+part*zion0(5,m)*(wpi(1,m)-vpara_dvi_a(m))       ! particle radial flux due to ExB
      ONE_me(16,j2)=ONE_me(16,j2)+(1.0_wp-part)*zion0(5,m)*(wpi(1,m)-vpara_dvi_a(m))

      ONE_me(17,j1)=ONE_me(17,j1)+part*zion0(5,m)*(wpi(1,m)-vpara_dvi_a(m))*energy ! energy radial flux due to ExB
      ONE_me(17,j2)=ONE_me(17,j2)+(1.0_wp-part)*zion0(5,m)*(wpi(1,m)-vpara_dvi_a(m))*energy

      vpara0=vpara         !  a bug fixed in April 2016

      !  calculate distribution functions at specified boxes
      !  NOTE:  need to work out OpenMP for this part 9/08 !!!
      do k=1,nprobe
         r_area=abs(zion0(1,m)-rprobe(k))
         if(r_area.le.dradius) then       ! within given radial range
            the1=abs(tprobe(k)-zion0(2,m))
            the2=abs(twopi+tprobe(k)-zion0(2,m))
            the_area=min(the1,the2)
            if(the_area.le.dtheta) then   ! within given theta range

               vper=sqrt(wpi(3,m)*2.0_wp)/vthprobe(k)  ! in local v_th unit
               vpara=vpara0
               vpara=vpara/vthprobe(k)            ! in local v_th unit
               if(vper.le.vmax.and.abs(vpara).le.vmax) then

                  yy=vpara/delta_v
                  mc=floor(yy)
                  k1=min(nv,max(-nv,mc))        !{-nv,nv}
                  k2=min(nv,max(-nv,mc+1))      !{-nv,nv}
                  !   part_k1=float(k2)-yy
                  part_k1=real(k2,wp)-yy

                  yy=vper/delta_v
                  mc=floor(yy)
                  l1=min(nv,max(0,mc))        !{0,nv}
                  l2=min(nv,max(0,mc+1))      !{0,nv}
                  !   part_l1=float(l2)-yy
                  part_l1=real(l2,wp)-yy

                  ! delta-f distriution
                  ss1=part_k1*part_l1*zion0(5,m)
                  ss2=part_k1*(1.0-part_l1)*zion0(5,m)
                  ss3=(1.0-part_k1)*part_l1*zion0(5,m)
                  ss4=(1.0-part_k1)*(1.0-part_l1)*zion0(5,m)

                  qq1=ss1
                  qq2=ss2
                  qq3=ss3
                  qq4=ss4

                  !                    df(k,k1,l1)=df(k,k1,l1)+part_k1*part_l1*zion0(5,m)
                  !                    df(k,k1,l2)=df(k,k1,l2)+part_k1*(1.0-part_l1)*zion0(5,m)
                  !                    df(k,k2,l1)=df(k,k2,l1)+(1.0-part_k1)*part_l1*zion0(5,m)
                  !                    df(k,k2,l2)=df(k,k2,l2)+(1.0-part_k1)*(1.0-part_l1)*zion0(5,m)

                  df(1,k,k1,l1)=df(1,k,k1,l1)+ss1
                  df(1,k,k1,l2)=df(1,k,k1,l2)+ss2
                  df(1,k,k2,l1)=df(1,k,k2,l1)+ss3
                  df(1,k,k2,l2)=df(1,k,k2,l2)+ss4

                  ! bootstrap current
                  df(7,k,k1,l1)=df(7,k,k1,l1)+ss1*vpara0*BB(m)
                  df(7,k,k1,l2)=df(7,k,k1,l2)+ss2*vpara0*BB(m)
                  df(7,k,k2,l1)=df(7,k,k2,l1)+ss3*vpara0*BB(m)
                  df(7,k,k2,l2)=df(7,k,k2,l2)+ss4*vpara0*BB(m)

                  ! particle flux distribution
                  df(6,k,k1,l1)=df(6,k,k1,l1)+ss1*wpi(1,m)
                  df(6,k,k1,l2)=df(6,k,k1,l2)+ss2*wpi(1,m)
                  df(6,k,k2,l1)=df(6,k,k2,l1)+ss3*wpi(1,m)
                  df(6,k,k2,l2)=df(6,k,k2,l2)+ss4*wpi(1,m)

                  ! momentum distribution
                  !                    am=zion0(4,m)*F(m)
                  ! replaced by:
                  am=Pa(m)*Pa(m)*dgf(m)
                  ss1=ss1*am
                  ss2=ss2*am
                  ss3=ss3*am
                  ss4=ss4*am

                  !                    df_m(k,k1,l1)=df_m(k,k1,l1)+part_k1*part_l1*zion0(5,m)*zion0(4,m)*F(m)
                  !                    df_m(k,k1,l2)=df_m(k,k1,l2)+part_k1*(1.0-part_l1)*zion0(5,m)*zion0(4,m)*F(m)
                  !                    df_m(k,k2,l1)=df_m(k,k2,l1)+(1.0-part_k1)*part_l1*zion0(5,m)*zion0(4,m)*F(m)
                  !                    df_m(k,k2,l2)=df_m(k,k2,l2)+(1.0-part_k1)*(1.0-part_l1)*zion0(5,m)*zion0(4,m)*F(m)

                  df(2,k,k1,l1)=df(2,k,k1,l1)+ss1
                  df(2,k,k1,l2)=df(2,k,k1,l2)+ss2
                  df(2,k,k2,l1)=df(2,k,k2,l1)+ss3
                  df(2,k,k2,l2)=df(2,k,k2,l2)+ss4

                  ! momentum flux distribution
                  ss1=ss1*wpi(1,m)
                  ss2=ss2*wpi(1,m)
                  ss3=ss3*wpi(1,m)
                  ss4=ss4*wpi(1,m)

                  !                    df_mf(k,k1,l1)=df_mf(k,k1,l1)+part_k1*part_l1*zion0(5,m)*zion0(4,m)*F(m)*wpi(1,m)
                  !                    df_mf(k,k1,l2)=df_mf(k,k1,l2)+part_k1*(1.0-part_l1)*zion0(5,m)*zion0(4,m)*F(m)*wpi(1,m)
                  !                    df_mf(k,k2,l1)=df_mf(k,k2,l1)+(1.0-part_k1)*part_l1*zion0(5,m)*zion0(4,m)*F(m)*wpi(1,m)
                  !                    df_mf(k,k2,l2)=df_mf(k,k2,l2)+(1.0-part_k1)*(1.0-part_l1)*zion0(5,m)*zion0(4,m)*F(m)*wpi(1,m)

                  df(3,k,k1,l1)=df(3,k,k1,l1)+ss1
                  df(3,k,k1,l2)=df(3,k,k1,l2)+ss2
                  df(3,k,k2,l1)=df(3,k,k2,l1)+ss3
                  df(3,k,k2,l2)=df(3,k,k2,l2)+ss4

                  ! energy distribution
                  qq1=qq1*energy
                  qq2=qq2*energy
                  qq3=qq3*energy
                  qq4=qq4*energy

                  df(4,k,k1,l1)=df(4,k,k1,l1)+qq1
                  df(4,k,k1,l2)=df(4,k,k1,l2)+qq2
                  df(4,k,k2,l1)=df(4,k,k2,l1)+qq3
                  df(4,k,k2,l2)=df(4,k,k2,l2)+qq4

                  ! energy flux distribution
                  qq1=qq1*wpi(1,m)
                  qq2=qq2*wpi(1,m)
                  qq3=qq3*wpi(1,m)
                  qq4=qq4*wpi(1,m)

                  df(5,k,k1,l1)=df(5,k,k1,l1)+qq1
                  df(5,k,k1,l2)=df(5,k,k1,l2)+qq2
                  df(5,k,k2,l1)=df(5,k,k2,l1)+qq3
                  df(5,k,k2,l2)=df(5,k,k2,l2)+qq4
               end if
            end if
         end if
      end do

   enddo
   !$omp end parallel

   call MPI_Reduce(ONE_me,TOTAL_pe,nn*nbin_a,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call MPI_Reduce(my_ion,total_ion,1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierr)

   call MPI_Reduce(df,df_tot,7*nprobe*(2*nv+1)*(nv+1),mpi_Rsize,MPI_SUM,1,MPI_COMM_WORLD,ierr)

   do i=1,nprobe
      df_heat(i,:,:)=df_tot(5,i,:,:)-5.0_wp/2.0_wp*Ti_box(i)*df_tot(6,i,:,:)  ! heat flux; added 2009
      df_tot(:,i,:,:)=df_tot(:,i,:,:)/sn_box(i)
      df_tot(:,i,:,:)=df_tot(:,i,:,:)*density_in_box(i)/(twopi*delta_v**2)   ! in unit of cm^(-3)*vth^(-3), not for all 5
      df_heat(i,:,:)=df_heat(i,:,:)/sn_box(i)
      df_heat(i,:,:)=df_heat(i,:,:)*density_in_box(i)/(twopi*delta_v**2)
      !     df_m_tot(i,:,:)=df_m_tot(i,:,:)/sn_box(i)
      !     df_m_tot(i,:,:)=df_m_tot(i,:,:)*density_in_box(i)/(twopi*delta_v**2)
      !     df_mf_tot(i,:,:)=df_mf_tot(i,:,:)/sn_box(i)
      !     df_mf_tot(i,:,:)=df_mf_tot(i,:,:)*density_in_box(i)/(twopi*delta_v**2)
   end do
   df_tot(7,:,:,:)=df_tot(7,:,:,:)*csc0_i   ! "bootstrap current"


   ! time=sum(TOTAL_pe(1,1:nbin_a))
   ! if(mype==0) print*, '+++++',(ONE_me(i,5),TOTAL_pe(i,5),i=1,4),time

   if(mype==0) then
      do kc=1,nbin_a

         !           wmean_0: delta_n/n_0
         !           dTi  : delta_T/T_0

         sum_n=1.0_wp/TOTAL_pe(1,kc)
         sum_sn0=1.0_wp/sn0(kc)

         !          wmean_0=TOTAL_pe(2,kc)*sum_n
         wmean_0=TOTAL_pe(2,kc)*sum_sn0  ! =delta_n/n_0
         w2mean_0=TOTAL_pe(4,kc)*sum_n

         av2=TOTAL_pe(3,kc)*sum_n
         Tg=2.0_wp*av2/factor*uveloci**2    ! in Kev

         delta_sumv2=2.0_wp*TOTAL_pe(5,kc)*sum_sn0

         flux_p=TOTAL_pe(6,kc)*sum_sn0*dnty_on_bin(kc)    ! in unit v_00, Gamma_r=int(dr/dt...)
         flux_E=TOTAL_pe(7,kc)*sum_sn0*dnty_on_bin(kc)    ! in unit E_00*v_00
         flux_vz=TOTAL_pe(8,kc)*sum_sn0
         flux_q=flux_E-5.0_wp/2.0_wp*Tr(kc)*flux_p          ! in unit E_00*v_00
         flux_PI=TOTAL_pe(10,kc)*sum_sn0*dnty_on_bin(kc)
         flux_PIp=TOTAL_pe(12,kc)*sum_sn0*dnty_on_bin(kc)   ! poloidal momentum flux
         flux_vp=TOTAL_pe(13,kc)*sum_sn0                    ! poloidal velocity
         angular_m=TOTAL_pe(9,kc)*sum_sn0*dnty_on_bin(kc)
         jb=TOTAL_pe(11,kc)*sum_sn0*dnty_on_bin(kc)*csc0_i  ! "boostrap current" in unit e0*v_00*B_00; 2010

         flux_p_strm=TOTAL_pe(14,kc)*sum_sn0*dnty_on_bin(kc)    ! in unit v_00, Gamma_r=int(dr/dt...)
         flux_E_strm=TOTAL_pe(15,kc)*sum_sn0*dnty_on_bin(kc)    ! in unit E_00*v_00
         flux_p_ExB =TOTAL_pe(16,kc)*sum_sn0*dnty_on_bin(kc)    ! in unit v_00, Gamma_r=int(dr/dt...)
         flux_E_ExB =TOTAL_pe(17,kc)*sum_sn0*dnty_on_bin(kc)    ! in unit E_00*v_00

         dTi=delta_sumv2/Tr(kc)/3.0_wp  ! =delta_T/T_0

         flux_p=flux_p*uveloci*1.0e4_wp         ! convert to 1/m^2.s
         flux_vz=flux_vz*uveloci*1.0e-5_wp      ! convert to km/s
         flux_vp=flux_vp*uveloci*1.0e-5_wp      ! convert to km/s
         flux_E=flux_E*uveloci*E_00*1.0e-3_wp   ! convert to J/m^2.s
         flux_q=flux_q*uveloci*E_00*1.0e-3_wp   ! convert to J/m^2.s

         flux_p_strm = flux_p_strm*uveloci*1.0e4_wp
         flux_E_strm = flux_E_strm*uveloci*E_00*1.0e-3_wp
         flux_p_ExB  = flux_p_ExB *uveloci*1.0e4_wp
         flux_E_ExB  = flux_E_ExB *uveloci*E_00*1.0e-3_wp

         write(nf(kc),800) istep,w2mean_0,wmean_0,TOTAL_pe(1,kc),Tg,dTi,flux_vz,angular_m, &
            flux_PI,flux_p,flux_E,flux_q,jb,total_ion/misum,flux_PIp,flux_vp, flux_p_strm, flux_E_strm, flux_p_ExB, flux_E_ExB
         call flush(nf(kc))

         w2_ion1(kc)=w2mean_0

      end do
   end if


   call mpi_bcast(w2_ion1,nbin_a,MPI_Rsize,0,MPI_COMM_WORLD,ierr)


   ONE_me_g=0.0_wp
   adum=0.0_wp
   !$omp parallel private(m,r,yy,kc,j1,j2,part,energy) reduction(+:ONE_me_g,adum)
   !$omp do
   do m=1,mi

      if (Flag_w_Threshold==1 .and. (zion0(5,m)<w_Thold_min .or. zion0(5,m)>w_Thold_max)) cycle

      r=zion0(7,m)
      yy=(r-s0)/da_size_g
      kc=int(yy)                    !{0,mpsi-1}
      kc=min(mpsi-1,max(0,kc))      !{0,mpsi-1}

      j1=kc
      j2=kc+1
      part=abs(real(j2,wp)-yy)

      vpara=zion0(4,m)*BB(m)
      energy=wpi(2,m)                ! v^2/2

      ONE_me_g(1,j1)=ONE_me_g(1,j1)+part*zion0(5,m)*energy         ! Temp contained in delta-f
      ONE_me_g(1,j2)=ONE_me_g(1,j2)+(1.0_wp-part)*zion0(5,m)*energy

      ONE_me_g(2,j1)=ONE_me_g(2,j1)+part*zion0(5,m)         ! density contained in delta-f
      ONE_me_g(2,j2)=ONE_me_g(2,j2)+(1.0_wp-part)*zion0(5,m)

      !        ONE_me_g(3,j1)=ONE_me_g(3,j1)+part*zion0(5,m)*wpi(1,m)*zion0(4,m)*F(m) ! angular momentum flux
      !        ONE_me_g(3,j2)=ONE_me_g(3,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*zion0(4,m)*F(m)
      ! replaced by:
      ONE_me_g(3,j1)=ONE_me_g(3,j1)+part*zion0(5,m)*wpi(1,m)*Pa(m)*Pa(m)*dgf(m) ! angular momentum flux
      ONE_me_g(3,j2)=ONE_me_g(3,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*Pa(m)*Pa(m)*dgf(m)

      ! added 2016
      ONE_me_g(10,j1)=ONE_me_g(10,j1)+part*zion0(5,m)*wpi(1,m)*Pa(m)*vpara ! parallel Reynolds stress associated to angular momentum flux
      ONE_me_g(10,j2)=ONE_me_g(10,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*Pa(m)*vpara

      !        ONE_me_g(4,j1)=ONE_me_g(4,j1)+part*zion0(5,m)*zion0(4,m)*F(m) ! angular momentum
      !        ONE_me_g(4,j2)=ONE_me_g(4,j2)+(1.0_wp-part)*zion0(5,m)*zion0(4,m)*F(m)
      ! replaced by:
      ONE_me_g(4,j1)=ONE_me_g(4,j1)+part*zion0(5,m)*Pa(m)*Pa(m)*dgf(m)  ! angular momentum
      ONE_me_g(4,j2)=ONE_me_g(4,j2)+(1.0_wp-part)*zion0(5,m)*Pa(m)*Pa(m)*dgf(m)


      ONE_me_g(5,j1)=ONE_me_g(5,j1)+part*zion0(5,m)*vpara ! v_// in delta-f
      ONE_me_g(5,j2)=ONE_me_g(5,j2)+(1.0_wp-part)*zion0(5,m)*vpara

      ONE_me_g(6,j1)=ONE_me_g(6,j1)+part*zion0(5,m)*vpara*BB(m) ! bootstrap current
      ONE_me_g(6,j2)=ONE_me_g(6,j2)+(1.0_wp-part)*zion0(5,m)*vpara*BB(m)

      ONE_me_g(7,j1)=ONE_me_g(7,j1)+part*zion0(5,m)*wpi(1,m)*TT(m)              ! poloidal angular momentum flux; modified Oct. 2012 with TT(m)= poloidal velocity
      ONE_me_g(7,j2)=ONE_me_g(7,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*TT(m)

      !        adum=adum+zion0(5,m)*zion0(4,m)*F(m)
      ! replaced by:
      adum=adum+zion0(5,m)*Pa(m)*Pa(m)*dgf(m)
      ! ZHIXIN: calculate particle and energy flux on grid
      ONE_me_g(8,j1)=ONE_me_g(8,j1)+part*zion0(5,m)*wpi(1,m) ! particle flux
      ONE_me_g(8,j2)=ONE_me_g(8,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)

      ONE_me_g(9,j1)=ONE_me_g(9,j1)+part*zion0(5,m)*wpi(1,m)*energy ! energy flux
      ONE_me_g(9,j2)=ONE_me_g(9,j2)+(1.0_wp-part)*zion0(5,m)*wpi(1,m)*energy

   enddo
   !$omp end parallel
   !ZHIXIN: size change from 7 to 9!!!
   call MPI_Reduce(ONE_me_g,TOTAL_pe_g,10*(mpsi+1),mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierr)
   call MPI_Reduce(adum,TM_sum,1,mpi_Rsize,MPI_SUM,numberpe-1,MPI_COMM_WORLD,ierr)
   if(mype==numberpe-1) then
      ONE_me_g(1,0:mpsi)=2.0_wp*TOTAL_pe_g(1,0:mpsi)/sn0_g(0:mpsi)/3.0_wp         ! =part of normalized <delta_T>

      TOTAL_pe_g(1,0:mpsi)=TOTAL_pe_g(1,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)
      TOTAL_pe_g(2,0:mpsi)=TOTAL_pe_g(2,0:mpsi)/sn0_g(0:mpsi)  ! =<delta_n>/n_0

      ONE_me_g(1,0:mpsi)=ONE_me_g(1,0:mpsi)-TOTAL_pe_g(2,0:mpsi)*vth_grid(0:mpsi)**2  ! =normalized <delta_T> -- corrected Aug. 2013

      TOTAL_pe_g(3,0:mpsi)=TOTAL_pe_g(3,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)    ! toroidal angular momentum flux
      TOTAL_pe_g(10,0:mpsi)=TOTAL_pe_g(10,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)  ! parallel Reynolds stress

      TOTAL_pe_g(4,0:mpsi)=TOTAL_pe_g(4,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)  ! angular momentum density
      ! in units of m_0*R_00*v_00*cm^(-3)
      TOTAL_pe_g(5,0:mpsi)=TOTAL_pe_g(5,0:mpsi)/sn0_g(0:mpsi)*vth_inv(0:mpsi)  ! (delta-v_//)/v_th
      TOTAL_pe_g(6,0:mpsi)=TOTAL_pe_g(6,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)*csc0_i  ! ion bootstrap current profile
      ! in unit e0*v_00*B_00; 2010
      TOTAL_pe_g(7,0:mpsi)=TOTAL_pe_g(7,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)  ! poloidal angular momentum flux

      enk_sum=0.0_wp
      do i=1,mpsi
         enk_sum=enk_sum+(TOTAL_pe_g(1,i-1)*jaco_sum(i-1)+TOTAL_pe_g(1,i)*jaco_sum(i))*  &
            (psi_grid(i)-psi_grid(i-1))*0.5_wp
      end do

      TOTAL_pe_g(1,0:mpsi)=ONE_me_g(1,0:mpsi)*vth_inv(0:mpsi)**2    !  <delta_T>/T_0 -- corrected Aug. 2013

      !ZHIXIN: particle and energy flux added
      TOTAL_pe_g(8,0:mpsi)=TOTAL_pe_g(8,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)  ! particle flux
      TOTAL_pe_g(9,0:mpsi)=TOTAL_pe_g(9,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)  ! energy flux

      TOTAL_pe_g(8,0:mpsi)=TOTAL_pe_g(8,0:mpsi)*uveloci*1.0e4_wp         ! convert to 1/m^2.s
      TOTAL_pe_g(9,0:mpsi)=TOTAL_pe_g(9,0:mpsi)*uveloci*E_00*1.0e-3_wp   ! convert to J/m^2.s

      write(93,802) istep, (TOTAL_pe_g(1,i),i=0,mpsi),enk_sum
      call flush(93)

      write(94,802) istep, (TOTAL_pe_g(2,i),i=0,mpsi)
      call flush(94)

      write(95,802) istep, (TOTAL_pe_g(3,i),i=0,mpsi)
      call flush(95)

      write(951,802) istep, (TOTAL_pe_g(10,i),i=0,mpsi)
      call flush(951)

      write(96,802) istep, (TOTAL_pe_g(4,i),i=0,mpsi),TM_sum
      call flush(96)

      write(97,802) istep, (TOTAL_pe_g(5,i),i=0,mpsi)
      call flush(97)

      write(99,802) istep, (TOTAL_pe_g(6,i),i=0,mpsi)
      call flush(99)

      write(85,802) istep, (TOTAL_pe_g(7,i),i=0,mpsi)
      call flush(85)
      !ZHIXIN add particle/energy flux
      write(86,802) istep, (TOTAL_pe_g(8,i),i=0,mpsi)
      call flush(86)

      write(87,802) istep, (TOTAL_pe_g(9,i),i=0,mpsi)
      call flush(87)

      ! direct numerical derivative to get d<delta_T>/da
      !$omp parallel do private(i,dpsi1,dpsi2,dpsi1_inv,dpsi2_inv,dpsi_inv,e_psi2,e_psi1)
      do i=1,mpsi-1
         dpsi1=psi_grid(i)-psi_grid(i-1)
         dpsi2=psi_grid(i+1)-psi_grid(i)
         dpsi1_inv=1.0_wp/dpsi1
         dpsi2_inv=1.0_wp/dpsi2
         dpsi_inv=1.0_wp/(dpsi1+dpsi2)
         e_psi2=(ONE_me_g(1,i+1)-ONE_me_g(1,i))*dpsi2_inv
         e_psi1=(ONE_me_g(1,i)-ONE_me_g(1,i-1))*dpsi1_inv
         ddT_da(i)=(e_psi1*dpsi2+e_psi2*dpsi1)*dpsi_inv
      enddo

      ! zero out d<delta_T>/da at boundary cells
      !    do i=0,nbound-1
      !       ddT_da(i)=ddT_da(i)*real(i,wp)/real(nbound,wp)
      !       ddT_da(mpsi-i)=ddT_da(mpsi-i)*real(i,wp)/real(nbound,wp)
      !    enddo
      do i=0,nbound
         rdum=smoothstep(0.0_wp, real(nbound,wp),real(i,wp))
         ddT_da(i)=ddT_da(i)*rdum
         ddT_da(mpsi-i)=ddT_da(mpsi-i)*rdum
      enddo

      dni_grid(:)=TOTAL_pe_g(2,:)
   end if

   call MPI_BCAST(ddT_da,mpsi+1,mpi_Rsize,numberpe-1,MPI_COMM_WORLD,ierr)

   call MPI_BCAST(dni_grid,mpsi+1,mpi_Rsize,numberpe-1,MPI_COMM_WORLD,ierr)

   if(mype.gt.0.and.mype.lt.(nfile+1)) then
      write(my_nf,801) istep, (phip00(i),i=my_ZL(1),my_ZL(2))
      call flush(my_nf)
      write(my_nfL,801) istep, (phi00(i)/Tc,i=my_ZL(1),my_ZL(2))
      call flush(my_nfL)
   end if

   if(mype==0) then
      write(80,802) istep, (fluxpsi(i),i=0,mpsi)
      call flush(80)

      write(83,802) istep, (flux_11(i),i=0,mpsi)
      call flush(83)
   end if

   factor=pmass0*uveloci**2*ulength/charg0**2/(4.0_wp*twopi)
   if(mype==numberpe-1) then
      En_tur=En_tur*factor  ! normalized field energy density, flux-surface-averaged
      En_tur_sum=0.0_wp
      do i=1,mpsi
         En_tur_sum=En_tur_sum+(En_tur(i-1)*jaco_sum(i-1)+En_tur(i)*jaco_sum(i))*  &
            (psi_grid(i)-psi_grid(i-1))*0.5_wp
      end do
      write(90,802) istep, (En_tur(i),i=0,mpsi),En_tur_sum
      call flush(90)
   end if

   if(mype==numberpe-2) then
      En_00=En_00*factor  ! normalized field energy density, flux-surface-averaged
      En_00_sum=0.0_wp
      do i=1,mpsi
         En_00_sum=En_00_sum+(En_00(i-1)*jaco_sum(i-1)+En_00(i)*jaco_sum(i))*  &
            (psi_grid(i)-psi_grid(i-1))*0.5_wp
      end do
      write(91,802) istep, (En_00(i),i=0,mpsi),En_00_sum
      call flush(91)
   end if

   if(mype==numberpe-4) then
      En_tot=En_tot*factor  ! normalized field energy density, flux-surface-averaged
      En_tot_sum=0.0_wp
      do i=1,mpsi
         En_tot_sum=En_tot_sum+(En_tot(i-1)*jaco_sum(i-1)+En_tot(i)*jaco_sum(i))*  &
            (psi_grid(i)-psi_grid(i-1))*0.5_wp
      end do
      write(92,802) istep, (En_tot(i),i=0,mpsi),En_tot_sum
      call flush(92)
   end if

   if(mype==numberpe-3) then
      write(98,802) istep, (ddT_da(i),i=0,mpsi)
      call flush(98)
   end if

   if(mype==1) then
      do k=1,nprobe
         write(201,400) istep, ( (df_tot(1,k,i,j), i=-nv,nv), j=0,nv)
         write(202,400) istep, ( (df_tot(2,k,i,j), i=-nv,nv), j=0,nv)
         write(203,400) istep, ( (df_tot(3,k,i,j), i=-nv,nv), j=0,nv)
         write(204,400) istep, ( (df_tot(4,k,i,j), i=-nv,nv), j=0,nv)
         write(205,400) istep, ( (df_tot(5,k,i,j), i=-nv,nv), j=0,nv)
         write(206,400) istep, ( (df_tot(6,k,i,j), i=-nv,nv), j=0,nv)
         write(207,400) istep, ( (df_heat(k,i,j),  i=-nv,nv), j=0,nv)
         write(208,400) istep, ( (df_tot(7,k,i,j), i=-nv,nv), j=0,nv)
      end do
      call flush_(201)
      call flush_(202)
      call flush_(203)
      call flush_(204)
      call flush_(205)
      call flush_(206)
      call flush_(207)
      call flush_(208)
   end if


   800 format(1x,i6,1x,30(d13.6,1x))
   801 format(1x,i6, 21(1x,d13.6))
   802 format(1x,i6, 503(1x,d13.6)) ! for maxium mpsi=500
   ! 802 format(1x,i6, 10(1x,d12.6)) ! for maxium mpsi=500
   400 format(1x,i6, 600(1x,d13.6))
   ! 400 format(1x,i6, 10(1x,d12.6))

end subroutine get_fluxes

