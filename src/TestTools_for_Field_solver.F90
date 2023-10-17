
module TestTools_for_Field_solver
   use precision, only:wp
   use Structured_Grid
   use MGY_PETSc_Interface_Data
   implicit none

   private

   type, public :: target_function

      integer :: m, n
      real(wp), allocatable :: sol_f_2D(:,:), src_2D(:,:), sim_f_2D(:,:) ! dimension=(SG%dim_gq, SG%dim_s)
      ! real :: test_src_3D(SG%dim_gq, SG%dim_s, mzetamax)
      contains
         procedure :: allocate_f_src_2D
         procedure :: deallocate_f_src_2D
         procedure :: set_mn_modes => construct_target_function
         procedure :: func_eval => evaluate_target_function
         procedure :: Jacob_Polar_T
         procedure :: Cal_local_highResolution_LHS_GK_Poisson
         procedure :: Cal_local_highResolution_LHS_GK_Poisson_parallelized
         procedure :: Cal_local_highResolution_LHS_GK_Poisson_parallelized_manual
         procedure :: test_GK_Poisson_solver
         procedure :: test_multiple_GK_Poisson_solvers
   end type target_function

   contains

      subroutine allocate_f_src_2D(this)
         implicit none
         class(target_function) this
         allocate(this%sol_f_2D(SG%dim_gq, SG%dim_s), this%src_2D(SG%dim_gq, SG%dim_s), this%sim_f_2D(SG%dim_gq, SG%dim_s))
      end subroutine allocate_f_src_2D

      subroutine deallocate_f_src_2D(this)
         implicit none
         class(target_function) this
         deallocate(this%sol_f_2D, this%src_2D, this%sim_f_2D)
      end subroutine deallocate_f_src_2D

      subroutine construct_target_function(this, in_m, in_n)
         class(target_function) :: this
         integer, intent(in), optional :: in_m,in_n

         this%m = 6
         this%n = 3

         if(present(in_m)) this%m = in_m
         if(present(in_n)) this%n = in_n

      end subroutine construct_target_function

      pure function evaluate_target_function(this,s,gq,gf,ddu_ij) result(out_val)
         implicit none

         class(target_function), intent(in) :: this
         integer, intent(in), optional :: ddu_ij ! derivative along i,j coordinate

         real(wp), intent(in) :: s,gq,gf
         real(wp) out_val


         if(present(ddu_ij)) then
            select case(ddu_ij)
               case(1) ! d/ds
                  out_val =  (2*pi/(s1-s0))*(sin(2*pi/(s1-s0)*(s-s0)))*cos(this%m*gq-this%n*gf)
               case(2) ! d/dgq 
                  out_val =  -this%m*(1.0-cos(2*pi/(s1-s0)*(s-s0)))*sin(this%m*gq-this%n*gf)
               case(3) ! d/dgf
                  out_val =  this%n*(1.0-cos(2*pi/(s1-s0)*(s-s0)))*sin(this%m*gq-this%n*gf)
               case(11) ! d^2/ds^2
                  out_val =  (2*pi/(s1-s0))**2.0*(cos(2*pi/(s1-s0)*(s-s0)))*cos(this%m*gq-this%n*gf)
               case(12, 21) ! d^2/(ds*dgq)
                  out_val = -this%m*(2*pi/(s1-s0))*(sin(2*pi/(s1-s0)*(s-s0)))*sin(this%m*gq-this%n*gf)
               case(13, 31) ! d^2/(ds*dgf)
                  out_val = this%n*(2*pi/(s1-s0))*(sin(2*pi/(s1-s0)*(s-s0)))*sin(this%m*gq-this%n*gf)
               case(22) ! d^2/dgq^2
                  out_val =  -this%m**2.0*(1.0-cos(2*pi/(s1-s0)*(s-s0)))*cos(this%m*gq-this%n*gf)
               case(23, 32) ! d^2/dgqdgf
                  out_val =  this%m*this%n*(1.0-cos(2*pi/(s1-s0)*(s-s0)))*cos(this%m*gq-this%n*gf)
               case(33) ! d^2/dgf^2
                  out_val =  -this%n**2.0*(1.0-cos(2*pi/(s1-s0)*(s-s0)))*cos(this%m*gq-this%n*gf)
            end select
         else
            out_val =  (1.0-cos(2*pi/(s1-s0)*(s-s0)))*cos(this%m*gq-this%n*gf)
         endif

         out_val = out_val/(uSI%phi0)
      end function evaluate_target_function


      pure subroutine Jacob_Polar_T(this, s,gq,gf, out_Jacob_Polar_T, ui, uj)
         use global_parameters, only: Flag_Def_dne_adiab
         use Magnetic_Grid
         implicit none

         class(target_function), intent(in) :: this
         real(wp), intent(in) :: s,gq,gf
         integer, optional, intent(in) :: ui, uj
         real(wp), intent(out) :: out_Jacob_Polar_T(3,3)
         real(wp) :: intp_a
         integer :: sid, gqid, i,j
         real(wp) :: ws0,ws1,wq0,wq1
         real(wp) :: ni_tot, ne_tot
         real(wp) :: intp_BB, invBsq, intp_B0_ctnv3(3)
         real(wp) :: loc_Metric_T(3,3)
         real(wp) :: intp_Jacob_sqf, tmp_factor
         real(wp) :: norm_gq, norm_gf

         norm_gq = modulo(gq,2.0*pi)
         norm_gf = modulo(gf,2.0*pi)

         sid = max(1,min(int((s-s0)*SG%inv_ds) +1,SG%dim_s-1)) ! +1 is due to Fortran's default index starting from 1
         ws1 = (s-SG%s_1D(sid))*SG%inv_ds
         ws0 = 1.0-ws1

         intp_a = ws0*SG%a_1D(sid) + ws1*SG%a_1D(sid+1)

         gqid = max(1,min(SG_dim_gq-1, int(norm_gq*SG_inv_dgq+1)))
         wq1 = (norm_gq-SG%gq_1D(gqid))*SG%inv_dgq
         wq0 = 1.0-wq1

         ni_tot = (ws0*SG%n0_1D(sid)+ws1*SG%n0_1D(sid+1))+ ws0*wq0*SG%dni_SM(gqid,sid)+ws0*wq1*SG%dni_SM(gqid+1,sid)+ws1*wq0*SG%dni_SM(gqid,sid+1)+ws1*wq1*SG%dni_SM(gqid+1,sid+1)
         ne_tot = (ws0*SG%n0_1D(sid)+ws1*SG%n0_1D(sid+1))+ ws0*wq0*SG%dne_SM(gqid,sid)+ws0*wq1*SG%dne_SM(gqid+1,sid)+ws1*wq0*SG%dne_SM(gqid,sid+1)+ws1*wq1*SG%dne_SM(gqid+1,sid+1)

         ! Density must be positive so ignore negative values because of numerical noises
         ne_tot = max(0.0,ne_tot)
         ni_tot = max(0.0,ni_tot)

         call TriCubic_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, BB, grad_BB, d23_BB, intp_BB)
         call TriCubic_Interpolation_of_Vec(intp_a,norm_gq,norm_gf, B_cntv, grad_B_cntv, d23_B_cntv, intp_B0_ctnv3)
         call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, Jacob_sqf, intp_Jacob_sqf)

         loc_Metric_T = 0.0
         out_Jacob_Polar_T = 0.0

         invBsq=1.0/(intp_BB**2.0) ! unit (T^2)
         tmp_factor=invBsq*(msm0*ne_tot+ msm0_i*ni_tot)

         if(present(ui) .and. present(uj)) then
            if(ui==1 .and. uj==1) then
               call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_ss_cntv, loc_Metric_T(1,1))
            elseif((ui==1 .and. uj==2) .or. (ui==2 .and. uj==1)) then
               call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_sq_cntv, loc_Metric_T(2,1))
               loc_Metric_T(1,2) = loc_Metric_T(2,1)
            elseif(ui==2 .and. uj==2) then
               call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_qq_cntv, loc_Metric_T(2,2))
            elseif(ui==3 .and. uj==3) then
               call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_ff_cntv, loc_Metric_T(3,3))
            endif

            out_Jacob_Polar_T(ui,uj) = tmp_factor*(loc_Metric_T(ui,uj)- invBsq*(intp_B0_ctnv3(ui)*intp_B0_ctnv3(uj)))
         else
            call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_ss_cntv, loc_Metric_T(1,1))
            call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_sq_cntv, loc_Metric_T(2,1))
            loc_Metric_T(1,2) = loc_Metric_T(2,1)
            call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_qq_cntv, loc_Metric_T(2,2))
            call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_ff_cntv, loc_Metric_T(3,3))

            do j=1,3
               do i=1,3
                  out_Jacob_Polar_T(i,j) = tmp_factor*(loc_Metric_T(i,j)- invBsq*(intp_B0_ctnv3(i)*intp_B0_ctnv3(j)))
               enddo
            enddo
         endif
         out_Jacob_Polar_T = intp_Jacob_sqf*out_Jacob_Polar_T

      end subroutine Jacob_Polar_T


      subroutine Cal_local_highResolution_LHS_GK_Poisson(this, out_LHS_f)
         implicit none
         class(target_function) :: this
         real(wp), intent(out) :: out_LHS_f(SG%dim_gq, SG%dim_s)
         real(wp) da,ds,dgq,dgf ! small value of delta_s, delta_gq, delta_gf
         real(wp) x0(3),dx_vec(3,3), pos_vec(3)
         integer :: ss,qq,ff
         integer :: i,j
         real(wp) Jacob
         real(wp) src
         real(wp) T1(3,3), T2(3,3)
         real(wp) f1,f2,f3,f4

         ds=1e-4
         dgq=1e-4
         dgf=0.001

         dx_vec(:,1)=[ds,  0.0, 0.0]
         dx_vec(:,2)=[0.0, dgq, 0.0]
         dx_vec(:,3)=[0.0, 0.0, dgf]

         out_LHS_f = 0.0

         do ss=2,SG%dim_s-1
            do qq=1,SG%dim_gq-1

               x0 = [SG%s_1D(ss),SG%gq_1D(qq),SG%zetamax]

               out_LHS_f(qq,ss) = 0.0
               do i=1,3
                  do j=1,3
                     pos_vec = x0 + 0.5*dx_vec(:,i)
                     call this%Jacob_Polar_T(pos_vec(1),pos_vec(2),pos_vec(3),T1)
                     pos_vec = x0 - 0.5*dx_vec(:,i)
                     call this%Jacob_Polar_T(pos_vec(1),pos_vec(2),pos_vec(3),T2)

                     pos_vec = x0 +0.5*dx_vec(:,i) +0.5*dx_vec(:,j)
                     f1 = this%func_eval(pos_vec(1),pos_vec(2),pos_vec(3))
                     pos_vec = x0 +0.5*dx_vec(:,i) -0.5*dx_vec(:,j)
                     f2 = this%func_eval(pos_vec(1),pos_vec(2),pos_vec(3))
                     pos_vec = x0 -0.5*dx_vec(:,i) +0.5*dx_vec(:,j)
                     f3 = this%func_eval(pos_vec(1),pos_vec(2),pos_vec(3))
                     pos_vec = x0 -0.5*dx_vec(:,i) -0.5*dx_vec(:,j)
                     f4 = this%func_eval(pos_vec(1),pos_vec(2),pos_vec(3))

                     out_LHS_f(qq,ss) = out_LHS_f(qq,ss) - 1.0/(SG%Jacob_sqf(qq,ss)*dx_vec(i,i)*dx_vec(j,j))*(T1(i,j)*(f1-f2)-T2(i,j)*(f3-f4))
                  enddo
               enddo
            enddo
         enddo

      end subroutine Cal_local_highResolution_LHS_GK_Poisson

      subroutine Cal_local_highResolution_LHS_GK_Poisson_parallelized(this, out_LHS_f)
         use particle_decomp
         use MGY_PETSc_Interface_Data
         implicit none
         class(target_function) :: this
         real(wp), intent(out) :: out_LHS_f(SG%dim_gq, SG%dim_s)
         real(wp) da,ds,dgq,dgf ! small value of delta_s, delta_gq, delta_gf
         real(wp) x0(3),dx_vec(3,3), pos_vec(3)
         integer :: ss,qq,ff
         integer :: i,j, Lid, Gid, ierr
         real(wp) Jacob
         real(wp) src
         real(wp) T1(3,3), T2(3,3)
         real(wp) f1,f2,f3,f4
         real(wp) LHS_f_mype(Lsize), LHS_f_partd(sum_Lsize_partd)

         ds=1e-3
         dgq=1e-3
         dgf=0.01

         dx_vec(:,1)=[ds,  0.0, 0.0]
         dx_vec(:,2)=[0.0, dgq, 0.0]
         dx_vec(:,3)=[0.0, 0.0, dgf]

         out_LHS_f = 0.0

         do Lid=1,Lsize
            Gid=mapL2G(Lid)
            qq=mapG2P(Gid)%gqid
            ss=mapG2P(Gid)%aid

            x0 = [SG%s_1D(ss),SG%gq_1D(qq),SG%zetamax]

            LHS_f_mype(Lid) = 0.0
            do i=1,3
               do j=1,3
                  pos_vec = x0 + 0.5*dx_vec(:,i)
                  call this%Jacob_Polar_T(pos_vec(1),pos_vec(2),pos_vec(3),T1)
                  pos_vec = x0 - 0.5*dx_vec(:,i)
                  call this%Jacob_Polar_T(pos_vec(1),pos_vec(2),pos_vec(3),T2)

                  pos_vec = x0 +0.5*dx_vec(:,i) +0.5*dx_vec(:,j)
                  f1 = this%func_eval(pos_vec(1),pos_vec(2),pos_vec(3))
                  pos_vec = x0 +0.5*dx_vec(:,i) -0.5*dx_vec(:,j)
                  f2 = this%func_eval(pos_vec(1),pos_vec(2),pos_vec(3))
                  pos_vec = x0 -0.5*dx_vec(:,i) +0.5*dx_vec(:,j)
                  f3 = this%func_eval(pos_vec(1),pos_vec(2),pos_vec(3))
                  pos_vec = x0 -0.5*dx_vec(:,i) -0.5*dx_vec(:,j)
                  f4 = this%func_eval(pos_vec(1),pos_vec(2),pos_vec(3))

                  LHS_f_mype(Lid) = LHS_f_mype(Lid) - 1.0/(SG%Jacob_sqf(qq,ss)*dx_vec(i,i)*dx_vec(j,j))*(T1(i,j)*(f1-f2)-T2(i,j)*(f3-f4))
               enddo
            enddo
         enddo

         call MPI_ALLGATHERV(LHS_f_mype,Lsize,mpi_Rsize,LHS_f_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
         do Lid=1,sum_Lsize_partd
            Gid=recv_mapL2G_partd(Lid)
            out_LHS_f(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = LHS_f_partd(Lid)
         enddo
         out_LHS_f(SG%dim_gq,:)=out_LHS_f(1,:)

      end subroutine Cal_local_highResolution_LHS_GK_Poisson_parallelized

      subroutine Cal_local_highResolution_LHS_GK_Poisson_parallelized_manual(this, out_LHS_f)
         use particle_decomp
         use MGY_PETSc_Interface_Data
         implicit none
         class(target_function) :: this
         real(wp), intent(out) :: out_LHS_f(SG%dim_gq, SG%dim_s)
         real(wp) ds,dgq,dgf ! small value of delta_s, delta_gq, delta_gf
         real(wp) x0(3)
         integer :: ss,qq,ff
         integer :: i,j, Lid, Gid, ierr
         real(wp) Jacob
         real(wp) src
         real(wp) T0(3,3), Tp(3,3), Tm(3,3)
         real(wp) LHS_f_mype(Lsize), LHS_f_partd(sum_Lsize_partd)

         ds=1e-3
         dgq=1e-3
         dgf=0.01

         out_LHS_f = 0.0

         do Lid=1,Lsize
            Gid=mapL2G(Lid)
            qq=mapG2P(Gid)%gqid
            ss=mapG2P(Gid)%aid

            x0 = [SG%s_1D(ss),SG%gq_1D(qq),SG%zetamax]

            LHS_f_mype(Lid) = 0.0
         
            call this%Jacob_Polar_T(x0(1), x0(2), x0(3), T0)

            ! J*eps^ss* d^2f/ds^2
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - T0(1,1)*this%func_eval(x0(1),x0(2),x0(3), 11)
            ! J*eps^gqgq * d^2f/dgq^2
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - T0(2,2)*this%func_eval(x0(1),x0(2),x0(3), 22)
            ! J*eps^gfgf * d^2f/dgf^2
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - T0(3,3)*this%func_eval(x0(1),x0(2),x0(3), 33)
            ! 2.0* J*eps^sgq* d^2f/dsdgq
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - 2.0*T0(1,2)*this%func_eval(x0(1),x0(2),x0(3), 12)
            ! 2.0* J*eps^sgf* d^2f/dsdgf
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - 2.0*T0(1,3)*this%func_eval(x0(1),x0(2),x0(3), 13)
            ! 2.0* J*eps^gqgf* d^2f/dgqdgf
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - 2.0*T0(2,3)*this%func_eval(x0(1),x0(2),x0(3), 23)


            call this%Jacob_Polar_T(x0(1)+0.5*ds, x0(2), x0(3), Tp)
            call this%Jacob_Polar_T(x0(1)-0.5*ds, x0(2), x0(3), Tm)
            ! d(J*eps^ss)/ds * df/ds
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - (Tp(1,1)-Tm(1,1))/ds*this%func_eval(x0(1),x0(2),x0(3), 1)
            ! d(J*eps^sgq)/ds * df/dgq
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - (Tp(1,2)-Tm(1,2))/ds*this%func_eval(x0(1),x0(2),x0(3), 2)
            ! d(J*eps^sgf)/ds * df/dgf
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - (Tp(1,3)-Tm(1,3))/ds*this%func_eval(x0(1),x0(2),x0(3), 3)


            call this%Jacob_Polar_T(x0(1), x0(2)+0.5*dgq, x0(3), Tp)
            call this%Jacob_Polar_T(x0(1), x0(2)-0.5*dgq, x0(3), Tm)
            ! d(J*eps^gqs)/dgq * df/ds
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - (Tp(2,1)-Tm(2,1))/dgq*this%func_eval(x0(1),x0(2),x0(3), 1)
            ! d(J*eps^gqgq)/dgq * df/dgq
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - (Tp(2,2)-Tm(2,2))/dgq*this%func_eval(x0(1),x0(2),x0(3), 2)
            ! d(J*eps^gqgf)/dgq * df/dgf
            LHS_f_mype(Lid) = LHS_f_mype(Lid) - (Tp(2,3)-Tm(2,3))/dgq*this%func_eval(x0(1),x0(2),x0(3), 3)


            LHS_f_mype(Lid) = LHS_f_mype(Lid)/SG%Jacob_sqf(qq,ss)
         enddo

         call MPI_ALLGATHERV(LHS_f_mype,Lsize,mpi_Rsize,LHS_f_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
         do Lid=1,sum_Lsize_partd
            Gid=recv_mapL2G_partd(Lid)
            out_LHS_f(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = LHS_f_partd(Lid)
         enddo
         out_LHS_f(SG%dim_gq,:)=out_LHS_f(1,:)

      end subroutine Cal_local_highResolution_LHS_GK_Poisson_parallelized_manual

      subroutine test_GK_Poisson_solver(this, output_fname)
         use mpi_defs
         use particle_decomp
         use hdf5
         use h5lt
         implicit none

         class(target_function) :: this
         character(len=*), optional, intent(in) :: output_fname
         real(wp) tt0,tt1
         integer ss,qq,ff
         integer m,n,nid,mid
         integer Nm,Nn
         real(wp) loc_err_sq, loc_sol_sq, all_err_sq, all_sol_sq
         real(wp),allocatable :: m_modes(:), n_modes(:), rel_err_2D(:,:)
         real(wp) in_gf
         character(128) :: filename
         integer :: h5err,file_exist
         integer(HSIZE_T) :: h5dims(5)
         integer(HID_T) file_id

         call this%allocate_f_src_2D


         m_modes = [(nid, nid=-200,200,20)]
         n_modes = [(nid, nid=-200,200,20)]

         Nm = size(m_modes)
         Nn = size(n_modes)

         allocate(rel_err_2D(Nm,Nn))


         do nid=1,Nn
            do mid=1,Nm

               m = m_modes(mid)
               n = n_modes(nid)

               call this%set_mn_modes(m,n)

               tt0=MPI_WTIME()
               do ss=1,SG%dim_s
                  do qq =1,SG%dim_gq
                     this%sol_f_2D(qq,ss) = this%func_eval(SG%s_1D(ss), SG%gq_1D(qq),SG%zetamax)
                  enddo
               enddo
               tt1=MPI_WTIME()
               if(mype==0) write(stdout,'(A,I0,A,ES11.3,A)') "mype=",mype,", In test_GK_Poisson_solver, Get sol_f_2D takes ",tt1-tt0," secs"

               tt0=MPI_WTIME()
               ! call this%Cal_local_highResolution_LHS_GK_Poisson(this%src_2D)
               call this%Cal_local_highResolution_LHS_GK_Poisson_parallelized(this%src_2D)
               ! call this%Cal_local_highResolution_LHS_GK_Poisson_parallelized_manual(this%src_2D)
                  
               tt1=MPI_WTIME()
               if(mype==0) write(stdout,'(A,I0,A,ES11.3,A)') "mype=",mype,", Cal_local_highResolution_LHS_GK_Poisson takes ",tt1-tt0," secs"

               call Solve_Poisson3D_on_SG(this%src_2D)

               this%sim_f_2D = SG%phi(:,:,1)

               loc_err_sq = sum((this%sim_f_2D-this%sol_f_2D)**2.0)
               loc_sol_sq = sum((this%sol_f_2D)**2.0)

               tt0=MPI_WTIME()
               call MPI_REDUCE(loc_err_sq, all_err_sq, 1, MPI_DOUBLE, MPI_SUM, 0, toroidal_comm,ierrr)
               call MPI_REDUCE(loc_sol_sq, all_sol_sq, 1, MPI_DOUBLE, MPI_SUM, 0, toroidal_comm,ierrr)
               tt1=MPI_WTIME()
               if(mype==0) write(stdout,'(A,I0,A,ES11.3,A)') "mype=",mype,", In test_GK_Poisson_solver, MPI_REDUCE takes ",tt1-tt0," secs"

               if(mype==0) then
                  write(stdout,'(A,I0,A,I0,A,3ES15.4)') "Rel. error of mode number (m,n)=(",m,", ", n,") = ", sqrt(all_err_sq), sqrt(all_sol_sq), sqrt(all_err_sq)/sqrt(all_sol_sq)
                  call flush(stdout)
               endif

               rel_err_2D(mid,nid) = sqrt(all_err_sq)/sqrt(all_sol_sq)
            enddo
         enddo

         do ss=1,SG%dim_s
            do qq=1,SG%dim_gq
               do ff=1,mzetamax
                  SG%test_data_3D(qq,ss,ff) = uSI%phi0*this%func_eval(SG%s_1D(ss), SG%gq_1D(qq), ff*(2.0)*pi/mzetamax)
               enddo
            enddo
         enddo


         if(mype==0) then
            ! Check if the file exists.
            filename = 'GK_Poisson_test_mn.h5'

            if(present(output_fname)) filename = output_fname

            inquire(file=filename, exist=file_exist)
            ! if (.not. file_exist) then
            !    write(stdout,'(A,A,A)'), "In test_GK_Poisson_solver, Error: File ", trim(filename), " does not exist."
            !    return
            ! endif

            call h5open_f(h5err)
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)

            h5dims(1) = Nm
            call h5ltmake_dataset_f(file_id, "m_modes",1,h5dims,H5T_NATIVE_DOUBLE, m_modes, h5err)
            h5dims(1) = Nn
            call h5ltmake_dataset_f(file_id, "n_modes",1,h5dims,H5T_NATIVE_DOUBLE, n_modes, h5err)

            h5dims(1) = SG%dim_s
            call h5ltmake_dataset_f(file_id, "skew_q",1,h5dims,H5T_NATIVE_DOUBLE, SG%skew_q_1D, h5err)

            h5dims(1:2) = [Nm,Nn]
            call h5ltmake_dataset_f(file_id, "rel_err_2D",2,h5dims,H5T_NATIVE_DOUBLE, rel_err_2D, h5err)

            call h5fclose_f(file_id, h5err)
            call h5close_f(h5err)
         endif
         call this%deallocate_f_src_2D

         deallocate(rel_err_2D)
      end subroutine test_GK_Poisson_solver


      subroutine test_multiple_GK_Poisson_solvers(this)
         implicit none

         class(target_function) :: this

         call Destruct_some_arrays_of_Field_Solver
         call Construct_GK_Poisson3D_Solver_using_Structured_Grid
         call this%test_GK_Poisson_solver("gk_poisson_mn_central.h5")

         call Destruct_some_arrays_of_Field_Solver
         call Construct_Metric_Tensor_on_locally_skewed_Structured_Grid
         call Construct_GK_Poisson3D_Solver_with_Skewed_Derivative
         call this%test_GK_Poisson_solver("gk_poisson_mn_skewed_q=natural.h5")

         call Destruct_some_arrays_of_Field_Solver
         SG%skew_q_1D = 2.0
         call Construct_Metric_Tensor_on_locally_skewed_Structured_Grid
         call Construct_GK_Poisson3D_Solver_with_Skewed_Derivative
         call this%test_GK_Poisson_solver("gk_poisson_mn_skewed_q=2.h5")

         call Destruct_some_arrays_of_Field_Solver
         SG%skew_q_1D = 1.0
         call Construct_Metric_Tensor_on_locally_skewed_Structured_Grid
         call Construct_GK_Poisson3D_Solver_with_Skewed_Derivative
         call this%test_GK_Poisson_solver("gk_poisson_mn_skewed_q=1.h5")

         call Destruct_some_arrays_of_Field_Solver
         SG%skew_q_1D = 0.5
         call Construct_Metric_Tensor_on_locally_skewed_Structured_Grid
         call Construct_GK_Poisson3D_Solver_with_Skewed_Derivative
         call this%test_GK_Poisson_solver("gk_poisson_mn_skewed_q=0.5.h5")


         call Destruct_some_arrays_of_Field_Solver
         SG%skew_q_1D = 0.2
         call Construct_Metric_Tensor_on_locally_skewed_Structured_Grid
         call Construct_GK_Poisson3D_Solver_with_Skewed_Derivative
         call this%test_GK_Poisson_solver("gk_poisson_mn_skewed_q=0.2.h5")

         call Destruct_some_arrays_of_Field_Solver
         call Construct_GK_Poisson3D_Solver_with_Extremal_Method
         call this%test_GK_Poisson_solver("gk_poisson_mn_extremal.h5")

      end subroutine test_multiple_GK_Poisson_solvers



end module TestTools_for_Field_solver