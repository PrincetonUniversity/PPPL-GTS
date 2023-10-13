
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
         procedure :: test_GK_Poisson_solver
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

      pure function evaluate_target_function(this,s,gq,gf) result(out_val)
         implicit none

         class(target_function) :: this

         real(wp), intent(in) :: s,gq,gf
         real(wp) out_val

         out_val =  (1.0-cos(2*pi/(s1-s0)*(s-s0)))*cos(this%m*gq-this%n*gf)
         out_val = out_val/(uSI%phi0)
      end function evaluate_target_function


      pure subroutine Jacob_Polar_T(this, s,gq,gf, out_Jacob_Polar_T)
         use global_parameters, only: Flag_Def_dne_adiab
         use Magnetic_Grid
         implicit none

         class(target_function) :: this
         real(wp), intent(in) :: s,gq,gf
         real(wp), intent(out) :: out_Jacob_Polar_T(3,3)
         real(wp) :: intp_a
         integer :: sid, gqid, i,j
         real(wp) :: ws0,ws1,wq0,wq1
         real(wp) :: ni_tot, ne_tot
         real(wp) :: intp_BB, invBsq, intp_B0_ctnv3(3)
         real(wp) :: loc_Metric_T(3,3)
         real(wp) :: intp_Jacob_sqf, tmp
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
         call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_ss_cntv, loc_Metric_T(1,1))
         call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_sq_cntv, loc_Metric_T(2,1))
         loc_Metric_T(1,2) = loc_Metric_T(2,1)
         call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_qq_cntv, loc_Metric_T(2,2))
         call TriLinear_Interpolation_of_Scalar(intp_a,norm_gq,norm_gf, g_ff_cntv, loc_Metric_T(3,3))

         out_Jacob_Polar_T = 0.0

         invBsq=1.0/(intp_BB**2.0) ! unit (T^2)
         do j=1,3
            do i=1,3
               tmp=invBsq*(msm0*ne_tot+ msm0_i*ni_tot)
               out_Jacob_Polar_T(i,j) = tmp*(loc_Metric_T(i,j)- invBsq*(intp_B0_ctnv3(i)*intp_B0_ctnv3(j)))
            enddo
         enddo

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

         ds=1e-3
         dgq=1e-3
         dgf=0.01

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


      subroutine test_GK_Poisson_solver(this)
         use mpi_defs
         use particle_decomp
         use hdf5
         use h5lt
         implicit none

         class(target_function) :: this
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

         call Destruct_some_arrays_of_Field_Solver
         ! call Construct_GK_Poisson3D_Solver_with_Skewed_Derivative
         call Construct_GK_Poisson3D_Solver_using_Structured_Grid

         call this%allocate_f_src_2D


         m_modes = [(nid, nid=-100,100,5)]
         n_modes = [(nid, nid=-100,100,5)]

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
               call this%Cal_local_highResolution_LHS_GK_Poisson(this%src_2D)
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

            h5dims(1:2) = [Nm,Nn]
            call h5ltmake_dataset_f(file_id, "rel_err_2D",2,h5dims,H5T_NATIVE_DOUBLE, rel_err_2D, h5err)

            call h5fclose_f(file_id, h5err)
            call h5close_f(h5err)
         endif
         call this%deallocate_f_src_2D

         deallocate(rel_err_2D)
      end subroutine test_GK_Poisson_solver

end module TestTools_for_Field_solver