subroutine get_fluxes_e(wpi)
  use precision
  use mpi_defs
  use global_parameters
  use particle_array
  use field_array
  use diagnosis_array
  use esi_interface
  
  implicit none

  real(wp)  wpi(3,me)
  integer   m, i, j, k, kc, j1, j2, m1, m2, l1, l2, mc, k1, k2,ierr
  real(wp)  adot,thedot,cosx,sinx,aj,da_table,yy,part,f_rho
  real(wp)  my_eon,total_eon

  real(wp), dimension((nn-2),nbin_a) :: ONE_me,TOTAL_pe
  real(wp), dimension(4,0:mpsi) :: ONE_me_g,TOTAL_pe_g

! add for diagnosis of distribution functions
  real(wp), dimension(7,nprobe,-nv:nv,0:nv) :: df, df_tot
                             ! df(1,:,:,:) --> distribution function
                             ! df(2,:,:,:) --> distribution of toroidal momentum
                             ! df(3,:,:,:) --> distribution of toroidal momentum flux
                             ! df(4,:,:,:) --> distribution of energy
                             ! df(5,:,:,:) --> distribution of energy flux
                             ! df(6,:,:,:) --> distribution of particle flux
                             ! df(7,:,:,:) --> distribution of "bootstrap current"
  real(wp), dimension(nprobe,-nv:nv,0:nv) :: df_heat   ! added 2009

  real(wp) r_area,vmax,the2,the1,part_l1,the_area,vper,ss1,ss2,ss3,ss4,am,qq1,qq2,qq3,qq4,upara0

  real(wp) factor,time,flux_p,flux_E,flux_q,flux_vz,flux_the,flux_PI,jb, &
         angular_m,diff_c,delta_sumv2,wmean_0,w2mean_0,w2mean_1, &
         vt,av2,Tg,dTi,avz,sum_n,sum_sn0,da_size,bootstrap,  &
         v2_lab, veloc, part_k1,  &
         En_tur_sum,En_00_sum,En_tot_sum,enk_sum,da_size_g

  real(wp) r,upara,energy

     factor=3000._wp*trans/pmass0
     f_rho=1.0_wp/for_rho   ! for_rho=(ms/m0)*(e0/es)

     da_size=(a1-a0)/real(nbin_a,wp)
     da_size_g=(s1-s0)/real(mpsi,wp)
     ONE_me=0.0_wp

     vmax=delta_v*real(nv)    ! maximum velocity in v_th unit
     df=0.0_wp
     df_tot=0.0_wp
     df_heat=0.0_wp
     my_eon=0.0_wp

!!     call esigetb(zeon0(1,:),zeon0(2,:),me)   ! get B
!     call esigetb(zeon0,me)   ! get B

!$omp parallel private(m,r,yy,kc,j1,j2,part,upara,upara0,energy) reduction(+:ONE_me,my_eon)
!$omp do
     do m=1,me
        r=zeon0(1,m)        
        yy=(r-a0)/da_size
        kc=int(yy+0.5_wp)            !{0,nbin_a}
        kc=min(nbin_a,max(0,kc))      !{0,nbin_a}

        j1=max(1,kc)
        j2=min(nbin_a,kc+1)
        part=abs(real(j2,wp)-0.5_wp-yy)

        upara=zeon0(4,m)*BB(m)*f_rho   ! normalized v_//
        energy=wpi(2,m)                ! v^2/2
     
        ONE_me(1,j1)=ONE_me(1,j1)+part
        ONE_me(1,j2)=ONE_me(1,j2)+(1.0_wp-part)

        ONE_me(2,j1)=ONE_me(2,j1)+part*zeon0(5,m)
        ONE_me(2,j2)=ONE_me(2,j2)+(1.0_wp-part)*zeon0(5,m)
       
        my_eon=my_eon+zeon0(5,m)

        ONE_me(3,j1)=ONE_me(3,j1)+part*energy
        ONE_me(3,j2)=ONE_me(3,j2)+(1.0_wp-part)*energy

        ONE_me(4,j1)=ONE_me(4,j1)+part*zeon0(5,m)**2
        ONE_me(4,j2)=ONE_me(4,j2)+(1.0_wp-part)*zeon0(5,m)**2

        ONE_me(5,j1)=ONE_me(5,j1)+part*zeon0(5,m)*energy         ! Temp contained in delta-f
        ONE_me(5,j2)=ONE_me(5,j2)+(1.0_wp-part)*zeon0(5,m)*energy

        ONE_me(6,j1)=ONE_me(6,j1)+part*zeon0(5,m)*wpi(1,m)       ! particle flux; wpi(1,m)=dr/dt|_EXB, not da/dt|_EXB
        ONE_me(6,j2)=ONE_me(6,j2)+(1.0_wp-part)*zeon0(5,m)*wpi(1,m)

        ONE_me(7,j1)=ONE_me(7,j1)+part*zeon0(5,m)*wpi(1,m)*energy ! energy flux
        ONE_me(7,j2)=ONE_me(7,j2)+(1.0_wp-part)*zeon0(5,m)*wpi(1,m)*energy

        ONE_me(8,j1)=ONE_me(8,j1)+part*zeon0(5,m)*upara ! v_// in delta-f
        ONE_me(8,j2)=ONE_me(8,j2)+(1.0_wp-part)*zeon0(5,m)*upara

        ONE_me(9,j1)=ONE_me(9,j1)+part*zeon0(5,m)*zeon0(4,m)*F(m)  ! angular momentum
        ONE_me(9,j2)=ONE_me(9,j2)+(1.0_wp-part)*zeon0(5,m)*zeon0(4,m)*F(m)

        ONE_me(10,j1)=ONE_me(10,j1)+part*zeon0(5,m)*wpi(1,m)*zeon0(4,m)*F(m)        ! angular momentum flux
        ONE_me(10,j2)=ONE_me(10,j2)+(1.0_wp-part)*zeon0(5,m)*wpi(1,m)*zeon0(4,m)*F(m)

        ONE_me(11,j1)=ONE_me(11,j1)+part*zeon0(5,m)*upara*BB(m) ! electron bootstrap current = v_//*B
        ONE_me(11,j2)=ONE_me(11,j2)+(1.0_wp-part)*zeon0(5,m)*upara*BB(m)

!  calculate distribution functions at specified boxes
!  NOTE:  need to work out OpenMP for this part 9/08 !!!
        do k=1,nprobe
           r_area=abs(zeon0(1,m)-rprobe(k))
           if(r_area.le.dradius) then       ! within given radial range
              the1=abs(tprobe(k)-zeon0(2,m))
              the2=abs(twopi+tprobe(k)-zeon0(2,m))
              the_area=min(the1,the2)
              if(the_area.le.dtheta) then   ! within given theta range

                 vper=sqrt(wpi(3,m)*m0ms*2.0_wp)/vthprobe_e(k)  ! in local v_th unit
                 upara0=upara
                 upara=upara/vthprobe_e(k)                      ! in local v_th unit
                 if(vper.le.vmax.and.abs(upara).le.vmax) then

                    yy=upara/delta_v
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
                    ss1=part_k1*part_l1*zeon0(5,m)
                    ss2=part_k1*(1.0-part_l1)*zeon0(5,m)
                    ss3=(1.0-part_k1)*part_l1*zeon0(5,m)
                    ss4=(1.0-part_k1)*(1.0-part_l1)*zeon0(5,m)

                    qq1=ss1
                    qq2=ss2
                    qq3=ss3
                    qq4=ss4

!                    df(k,k1,l1)=df(k,k1,l1)+part_k1*part_l1*zeon0(5,m)
!                    df(k,k1,l2)=df(k,k1,l2)+part_k1*(1.0-part_l1)*zeon0(5,m)
!                    df(k,k2,l1)=df(k,k2,l1)+(1.0-part_k1)*part_l1*zeon0(5,m)
!                    df(k,k2,l2)=df(k,k2,l2)+(1.0-part_k1)*(1.0-part_l1)*zeon0(5,m)

                    df(1,k,k1,l1)=df(1,k,k1,l1)+ss1
                    df(1,k,k1,l2)=df(1,k,k1,l2)+ss2
                    df(1,k,k2,l1)=df(1,k,k2,l1)+ss3
                    df(1,k,k2,l2)=df(1,k,k2,l2)+ss4

! bootstrap current
                    df(7,k,k1,l1)=df(7,k,k1,l1)+ss1*upara0*BB(m)
                    df(7,k,k1,l2)=df(7,k,k1,l2)+ss2*upara0*BB(m)
                    df(7,k,k2,l1)=df(7,k,k2,l1)+ss3*upara0*BB(m)
                    df(7,k,k2,l2)=df(7,k,k2,l2)+ss4*upara0*BB(m)

! particle flux distribution
                    df(6,k,k1,l1)=df(6,k,k1,l1)+ss1*wpi(1,m)
                    df(6,k,k1,l2)=df(6,k,k1,l2)+ss2*wpi(1,m)
                    df(6,k,k2,l1)=df(6,k,k2,l1)+ss3*wpi(1,m)
                    df(6,k,k2,l2)=df(6,k,k2,l2)+ss4*wpi(1,m)

! momentum distribution
                    am=zeon0(4,m)*F(m)
                    ss1=ss1*am
                    ss2=ss2*am
                    ss3=ss3*am
                    ss4=ss4*am

!                    df_m(k,k1,l1)=df_m(k,k1,l1)+part_k1*part_l1*zeon0(5,m)*zeon0(4,m)*F(m)
!                    df_m(k,k1,l2)=df_m(k,k1,l2)+part_k1*(1.0-part_l1)*zeon0(5,m)*zeon0(4,m)*F(m)
!                    df_m(k,k2,l1)=df_m(k,k2,l1)+(1.0-part_k1)*part_l1*zeon0(5,m)*zeon0(4,m)*F(m)
!                    df_m(k,k2,l2)=df_m(k,k2,l2)+(1.0-part_k1)*(1.0-part_l1)*zeon0(5,m)*zeon0(4,m)*F(m)

                    df(2,k,k1,l1)=df(2,k,k1,l1)+ss1
                    df(2,k,k1,l2)=df(2,k,k1,l2)+ss2
                    df(2,k,k2,l1)=df(2,k,k2,l1)+ss3
                    df(2,k,k2,l2)=df(2,k,k2,l2)+ss4

! momentum flux distribution
                    ss1=ss1*wpi(1,m)
                    ss2=ss2*wpi(1,m)
                    ss3=ss3*wpi(1,m)
                    ss4=ss4*wpi(1,m)

!                    df_mf(k,k1,l1)=df_mf(k,k1,l1)+part_k1*part_l1*zeon0(5,m)*zeon0(4,m)*F(m)*wpi(1,m)
!                    df_mf(k,k1,l2)=df_mf(k,k1,l2)+part_k1*(1.0-part_l1)*zeon0(5,m)*zeon0(4,m)*F(m)*wpi(1,m)
!                    df_mf(k,k2,l1)=df_mf(k,k2,l1)+(1.0-part_k1)*part_l1*zeon0(5,m)*zeon0(4,m)*F(m)*wpi(1,m)
!                    df_mf(k,k2,l2)=df_mf(k,k2,l2)+(1.0-part_k1)*(1.0-part_l1)*zeon0(5,m)*zeon0(4,m)*F(m)*wpi(1,m)

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

     ONE_me(9,:)=ONE_me(9,:)/for_rho      ! using normalized v_para/B <-- rho_para in momentum momentum flux 
     ONE_me(10,:)=ONE_me(10,:)/for_rho
  
  call MPI_Reduce(ONE_me,TOTAL_pe,(nn-2)*nbin_a,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(my_eon,total_eon,1,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  call MPI_Reduce(df,df_tot,7*nprobe*(2*nv+1)*(nv+1),mpi_Rsize,MPI_SUM,1,MPI_COMM_WORLD,ierr)

  do i=1,nprobe
     df_heat(i,:,:)=df_tot(5,i,:,:)-5.0_wp/2.0_wp*Te_box(i)*df_tot(6,i,:,:)  ! heat flux; added 2009
     df_tot(:,i,:,:)=df_tot(:,i,:,:)/sn_box(i)
     df_tot(:,i,:,:)=df_tot(:,i,:,:)*density_in_box(i)/(twopi*delta_v**2)   ! in unit of cm^(-3)*vth^(-3), not for all 5
     df_heat(i,:,:)=df_heat(i,:,:)/sn_box(i)
     df_heat(i,:,:)=df_heat(i,:,:)*density_in_box(i)/(twopi*delta_v**2)
!     df_m_tot(i,:,:)=df_m_tot(i,:,:)/sn_box(i)
!     df_m_tot(i,:,:)=df_m_tot(i,:,:)*density_in_box(i)/(twopi*delta_v**2)
!     df_mf_tot(i,:,:)=df_mf_tot(i,:,:)/sn_box(i)
!     df_mf_tot(i,:,:)=df_mf_tot(i,:,:)*density_in_box(i)/(twopi*delta_v**2)
  end do
     df_tot(7,:,:,:)=df_tot(7,:,:,:)*csc0   ! "bootstrap current"

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
          flux_q=flux_E-5.0_wp/2.0_wp*Tr_e(kc)*flux_p          ! in unit E_00*v_00
          flux_PI=TOTAL_pe(10,kc)*sum_sn0*dnty_on_bin(kc)
          angular_m=TOTAL_pe(9,kc)*sum_sn0*dnty_on_bin(kc)
          jb=TOTAL_pe(11,kc)*sum_sn0*dnty_on_bin(kc)*csc0  ! "boostrap current" in unit e0*v_00*B_00; 2010

          dTi=delta_sumv2/Tr_e(kc)/3.0_wp  ! =delta_T/T_0

          flux_p=flux_p*uveloci*1.0e4_wp         ! convert to 1/m^2.s
          flux_vz=flux_vz*uveloci*1.0e-5_wp      ! convert to km/s
          flux_E=flux_E*uveloci*E_00*1.0e-3_wp   ! convert to J/m^2.s
          flux_q=flux_q*uveloci*E_00*1.0e-3_wp   ! convert to J/m^2.s

          write(nf_e(kc),800) istep,w2mean_0,wmean_0,TOTAL_pe(1,kc),Tg,dTi,flux_vz,angular_m, &
                                  flux_PI,flux_p,flux_E,flux_q,jb,total_eon/mesum
          call flush(nf_e(kc))

          w2_eon1(kc)=w2mean_0

       end do
    end if


    call mpi_bcast(w2_eon1,nbin_a,MPI_Rsize,0,MPI_COMM_WORLD,ierr)


     ONE_me_g=0.0_wp
!$omp parallel private(m,r,yy,kc,j1,j2,part,energy) reduction(+:ONE_me_g)
!$omp do
     do m=1,me
        r=zeon0(7,m)
        yy=(r-s0)/da_size_g
        kc=int(yy)                  !{0,mpsi-1}
        kc=min(mpsi-1,max(0,kc))      !{0,mpsi-1}

        j1=kc
        j2=kc+1
        part=abs(real(j2,wp)-yy)

        upara=zeon0(4,m)*BB(m)*f_rho   ! normalized v_//
        energy=wpi(2,m)                ! v^2/2

        ONE_me_g(1,j1)=ONE_me_g(1,j1)+part*zeon0(5,m)*energy         ! Temp contained in delta-f
        ONE_me_g(1,j2)=ONE_me_g(1,j2)+(1.0_wp-part)*zeon0(5,m)*energy

        ONE_me_g(2,j1)=ONE_me_g(2,j1)+part*zeon0(5,m)         ! density contained in delta-f
        ONE_me_g(2,j2)=ONE_me_g(2,j2)+(1.0_wp-part)*zeon0(5,m)

        ONE_me_g(3,j1)=ONE_me_g(3,j1)+part*zeon0(5,m)*upara*BB(m) ! bootstrap current
        ONE_me_g(3,j2)=ONE_me_g(3,j2)+(1.0_wp-part)*zeon0(5,m)*upara*BB(m)

        ONE_me_g(4,j1)=ONE_me_g(4,j1)+part*zeon0(5,m)*wpi(1,m)*zeon0(4,m)*F(m) ! angular momentum flux
        ONE_me_g(4,j2)=ONE_me_g(4,j2)+(1.0_wp-part)*zeon0(5,m)*wpi(1,m)*zeon0(4,m)*F(m)

     enddo
!$omp end parallel

  call MPI_Reduce(ONE_me_g,TOTAL_pe_g,4*(mpsi+1),MPI_REAL8,MPI_SUM,numberpe-3,MPI_COMM_WORLD,ierr)
  if(mype==numberpe-3) then 
     ONE_me_g(1,0:mpsi)=2.0_wp*TOTAL_pe_g(1,0:mpsi)/sn0_g(0:mpsi)/3.0_wp         ! =part of normalized <delta_Te>

     TOTAL_pe_g(1,0:mpsi)=TOTAL_pe_g(1,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)
     TOTAL_pe_g(2,0:mpsi)=TOTAL_pe_g(2,0:mpsi)/sn0_g(0:mpsi)  ! =<delta_n>/n_0

     TOTAL_pe_g(3,0:mpsi)=TOTAL_pe_g(3,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)*csc0    ! electron bootstrap current profile
     TOTAL_pe_g(4,0:mpsi)=TOTAL_pe_g(4,0:mpsi)*dnty_on_bin_g(0:mpsi)/sn0_g(0:mpsi)/for_rho ! angular momentum flux

     enk_sum=0.0_wp
     do i=1,mpsi
        enk_sum=enk_sum+(TOTAL_pe_g(1,i-1)*jaco_sum(i-1)+TOTAL_pe_g(1,i)*jaco_sum(i))*  &
                (psi_grid(i)-psi_grid(i-1))*0.5_wp
     end do

     TOTAL_pe_g(1,0:mpsi)=ONE_me_g(1,0:mpsi)/Te_grid(0:mpsi)-TOTAL_pe_g(2,0:mpsi)  ! = <delta_Te>/Te_0

     write(193,802) istep, (TOTAL_pe_g(1,i),i=0,mpsi),enk_sum
     call flush(193)

     write(194,802) istep, (TOTAL_pe_g(2,i),i=0,mpsi)
     call flush(194)

     write(195,802) istep, (TOTAL_pe_g(4,i),i=0,mpsi)
     call flush(195)

     write(199,802) istep, (TOTAL_pe_g(3,i),i=0,mpsi)
     call flush(199)
  end if

    if(mype==1) then
       do k=1,nprobe
          write(301,400) istep, ( (df_tot(1,k,i,j), i=-nv,nv), j=0,nv)
          write(302,400) istep, ( (df_tot(2,k,i,j), i=-nv,nv), j=0,nv)
          write(303,400) istep, ( (df_tot(3,k,i,j), i=-nv,nv), j=0,nv)
          write(304,400) istep, ( (df_tot(4,k,i,j), i=-nv,nv), j=0,nv)
          write(305,400) istep, ( (df_tot(5,k,i,j), i=-nv,nv), j=0,nv)
          write(306,400) istep, ( (df_tot(6,k,i,j), i=-nv,nv), j=0,nv)
          write(307,400) istep, ( (df_heat(k,i,j),  i=-nv,nv), j=0,nv)
          write(308,400) istep, ( (df_tot(7,k,i,j), i=-nv,nv), j=0,nv)
       end do
       call flush_(301)
       call flush_(302)
       call flush_(303)
       call flush_(304)
       call flush_(305)
       call flush_(306)
       call flush_(307)
       call flush_(308)
    end if


800 format(1x,i6,1x,16(d12.6,1x))
801 format(1x,i6, 21(1x,d12.6))
802 format(1x,i6, 503(1x,d12.6)) ! for maxium mpsi=500
! 802 format(1x,i6, 10(1x,d12.6)) ! for maxium mpsi=500
400 format(1x,i6, 600(1x,d12.6))
! 400 format(1x,i6, 10(1x,d12.6))

end subroutine get_fluxes_e

