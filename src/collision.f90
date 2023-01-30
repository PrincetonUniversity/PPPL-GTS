subroutine collision
   use mpi_defs
   use global_parameters
   use particle_array
   use field_array
   use esi_interface
   use function_gg_module

   implicit none

   integer i

   !    NOTE: output from call n_T_u_group is omega_t and d(omega_t)/da;  12/07
   !          B(=BB) and I(=F) at particles are needed for calculating u_//=I*omega_t/B;
   !          BB, I can be obtained after call gcmotion, or call esigetbf, esiget2dfunctions.

   call n_T_u_group  ! prepare n, T, dlnn/da, dlnT/da, etc., which are also used in PUSHI/PUSHION

   if(icollision.gt.0) then   ! include collision
      i=mod(istep,icollision)
      if(i==0) then           ! calculate collisions every icollision-step
         if(irk==1) then      ! calculate collisions at irk=1

            ! Update needed particle information from DEVICE to HOST
            !$acc update self(backn(:mi),backu(:mi),backT(:mi))
            if(acc_newer_D_zion) then
               !$acc update self(zion(:,:mi))
            endif

            !           call esigetb(zion(1,:),zion(2,:),mi)   ! get B
            !           call esigetb(zion,mi)   ! get B-->BB
            call esigetbf(zion,mi,ziondim)   ! get B-->BB, I-->F  12/07

            !  output backu from call n_T_u_group is omega_t
            backu(1:mi)=backu(1:mi)*F(1:mi)/BB(1:mi) ! =u_//=I*omega_t/B

            call CTPandCELL

            !  make sure output backu=omega_t after exit this subroutine
            backu(1:mi)=backu(1:mi)*BB(1:mi)/F(1:mi) ! =omega_t


            !$acc update device(zion(:,:mi),backu(:mi))
            call n_T_u_group ! Update back_f0_i

          end if
      end if
   end if

end subroutine collision
