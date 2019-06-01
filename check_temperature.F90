subroutine check_temperature
   use precision
   use mpi_defs
   use global_parameters
   use particle_array
   use field_array
   use diagnosis_array
   use esi_interface

   implicit none

   integer,parameter :: mbin_psi=nbin_a
   integer gcmotion, m, i, j, k, kc, j1, j2, m1, m2, l1, l2, mc, k1, k2,ierr
   real(wp)  adot,thedot,cosx,sinx,aj,da_table,yy,part

   integer, dimension(mbin_psi) :: mf, mfL
   real(wp),  dimension(4,mbin_psi) :: ONE_me,TOTAL_pe

   real(wp) factor,time,flux_p,flux_E,flux_q,flux_vz,flux_the,flux_PI, &
      angular_m,diff_c,delta_sumv2,wmean_0,w2mean_0,w2mean_1, &
      vt,av2,Tg,dTi,avz,sum_n,sum_sn0,da_size,bootstrap,  &
      v2_lab, veloc, partv_j1, partv_j2, part_k1

   real(wp), dimension(10) :: aaa,the,Fj,Faj,gFaj,gYaj,rj,raj,rqj,zj,zaj,zqj,Bj

   real(wp) r,upara,energy

   do i=1,mbin_psi
      mf(i)=700+i
   end do

   factor=3000._wp*trans/pmass0

   da_size=(a1-a0)/real(mbin_psi,wp)
   ONE_me=0._wp

   !     call esigetb(zeon(1,:),zeon(2,:),mi)   ! get B
   call esigetb(zeon,mi,zeondim)   ! get B

   do m=1,mi
      r=zeon(1,m)
      yy=(r-a0)/da_size
      kc=int(yy+0.5_wp)            !{0,mbin_psi}
      kc=min(mbin_psi,max(0,kc))   !{0,mbin_psi}

      j1=max(1,kc)
      j2=min(mbin_psi,kc+1)
      part=abs(real(j2,wp)-0.5_wp-yy)

      upara=zeon(4,m)*BB(m)
      energy=max(1.0e-20_wp,upara*upara+zeon(6,m)*zeon(6,m)*BB(m)*2.0_wp)

      ONE_me(1,j1)=ONE_me(1,j1)+part
      ONE_me(1,j2)=ONE_me(1,j2)+(1.0_wp-part)

      ONE_me(2,j1)=ONE_me(2,j1)+part*zeon(5,m)
      ONE_me(2,j2)=ONE_me(2,j2)+(1.0_wp-part)*zeon(5,m)

      ONE_me(3,j1)=ONE_me(3,j1)+part*energy
      ONE_me(3,j2)=ONE_me(3,j2)+(1.0_wp-part)*energy

      ONE_me(4,j1)=ONE_me(4,j1)+part*zeon(5,m)**2
      ONE_me(4,j2)=ONE_me(4,j2)+(1.0_wp-part)*zeon(5,m)**2
   enddo


   call MPI_Reduce(ONE_me,TOTAL_pe,4*mbin_psi,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierr)


   ! time=sum(TOTAL_pe(1,1:mbin_psi))
   ! if(mype==0) print*, '+++++',(ONE_me(i,5),TOTAL_pe(i,5),i=1,4),time

   if(mype==0) then
      do kc=1,mbin_psi

         sum_n=1.0_wp/TOTAL_pe(1,kc)
         sum_sn0=1.0_wp/sn0(kc)

         wmean_0=TOTAL_pe(2,kc)*sum_n
         w2mean_0=TOTAL_pe(4,kc)*sum_n

         av2=TOTAL_pe(3,kc)*sum_sn0
         Tg=av2/factor*uveloci**2    ! in Kev

         TOTAL_pe(1,kc)=TOTAL_pe(1,kc)*sum_sn0

         write(mf(kc),800) istep,w2mean_0,wmean_0,TOTAL_pe(1,kc),Tg
         call flush(mf(kc))

      end do
   end if

   800 format(1x,i4,1x,10(ES12.6,1x))


end subroutine check_temperature

