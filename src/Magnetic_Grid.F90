module Magnetic_Grid
   use precision
   ! use global_parameters
   implicit none

   type esi_data
      integer dim_a,dim_gq !number of grid points in radial and poloidal direction
      real(wp),allocatable,dimension(:) :: a_1D, gq_1D
      real(wp),allocatable,dimension(:) :: F,Fa,gFa,gFaa,gYa,gYaa
      real(wp),allocatable,dimension(:) :: T,Ta,P,Pa
      real(wp),allocatable,dimension(:,:) :: R,Ra,Rq,Raq
      real(wp),allocatable,dimension(:,:) :: Z,Za,Zq,Zaq
      real(wp),allocatable,dimension(:,:) :: BB,BBa,BBq,BBaq
      real(wp),allocatable,dimension(:,:) :: gh,gh_a,gh_q,gh_aq
      real(wp),allocatable,dimension(:) :: ne,dneda

      !post-processed varaiables
      real(wp) R0,Z0
      real(wp),allocatable,dimension(:) :: gY,gY_norm,gF,gF_norm,safety_q

      real(wp),allocatable,dimension(:,:) :: Jacobian

   end type esi_data

   type FieldMap_3D
      integer dim_a,dim_gq,dim_gf
      real(wp) inv_da, inv_dgq, inv_dgf
      real(wp),allocatable,dimension(:) :: a_1D, gq_1D, gf_1D
      real(wp),allocatable,dimension(:,:,:) :: Max_BB0_CW, Max_BB0_CCW,BB0,Mirror_Ratio
   end type FieldMap_3D


   type(FieldMap_3D), target :: FMap3D ! This will be loaded by reading some pre-calculated hdf5 file
   type(esi_data), target :: esi
   integer dim_a,dim_gq,dim_gf ! Number of grid points
   !$acc declare create(dim_a,dim_gq,dim_gf)
   real(wp) da,dgq,dgf,inv_da,inv_dgq,inv_dgf
   !$acc declare create(da,dgq,dgf,inv_da,inv_dgq,inv_dgf)
   real(wp),allocatable,dimension(:) :: a_1D, gq_1D, gf_1D
   !$acc declare create(a_1D,gq_1D,gf_1D)
   real(wp),allocatable,dimension(:,:) :: R_2D, Z_2D  ! (R,Z) based on 2D equilibrium
   real(wp),allocatable,dimension(:,:,:) :: d12_R_2D, d12_Z_2D ! (R,Z)'s 1st and 2nd derivatvies in (a,gq) direction
   real(wp),allocatable,dimension(:,:,:) :: Jacobian_3D

   real(wp),allocatable,dimension(:,:,:) :: g_aa_cov,g_aq_cov,g_qq_cov,g_ff_cov



   ! Equilibrium B fields
   real(wp),allocatable,dimension(:,:,:) :: BB0 ! magnetic of B fields
   real(wp),allocatable,dimension(:,:,:,:) :: grad_BB0, d23_BB0
   real(wp),allocatable,dimension(:,:,:,:) :: B0_cov, B0_cntv, j0_cov, j0_cntv
   real(wp),allocatable,dimension(:,:,:,:,:) :: grad_B0_cov, grad_B0_cntv, grad_j0_cov, grad_j0_cntv
   real(wp),allocatable,dimension(:,:,:,:,:) :: d23_B0_cov, d23_B0_cntv, d23_j0_cov, d23_j0_cntv ! Higher order derivatvies

   ! dB (pertrubed) B fields
   real(wp),allocatable,dimension(:,:,:) :: dBB
   real(wp),allocatable,dimension(:,:,:,:) :: grad_dBB, d23_dBB
   real(wp),allocatable,dimension(:,:,:,:) :: dB_cov, dB_cntv, dj_cov, dj_cntv
   real(wp),allocatable,dimension(:,:,:,:,:) :: grad_dB_cov, grad_dB_cntv, grad_dj_cov, grad_dj_cntv
   real(wp),allocatable,dimension(:,:,:,:,:) :: d23_dB_cov, d23_dB_cntv, d23_dj_cov, d23_dj_cntv

   ! total B(=B0+dB)
   real(wp),allocatable,dimension(:,:,:) :: BB
   real(wp),allocatable,dimension(:,:,:,:) :: grad_BB, d23_BB
   real(wp),allocatable,dimension(:,:,:,:) :: B_cov,B_cntv,j_cov,j_cntv
   real(wp),allocatable,dimension(:,:,:,:,:) :: grad_B_cov, grad_B_cntv, grad_j_cov, grad_j_cntv
   real(wp),allocatable,dimension(:,:,:,:,:) :: d23_B_cov, d23_B_cntv, d23_j_cov, d23_j_cntv

   !$acc declare create(R_2D,d12_R_2D,Z_2D,d12_Z_2D)
   !$acc declare create(BB,grad_BB,d23_BB)
   !$acc declare create(B_cov,grad_B_cov,d23_B_cov)
   !$acc declare create(B_cntv,grad_B_cntv,d23_B_cntv)
   !$acc declare create(j_cntv,grad_j_cntv,d23_j_cntv)

contains

   subroutine Read_FMap3D_on_MG
      use hdf5
      use h5lt
      use global_parameters, only:stdout,mype
      implicit none
      type(FieldMap_3D), pointer :: FM
      character(128) :: filename
      integer :: Na,Nq,Nf,i,j,k
      integer :: h5err
      integer(HSIZE_T) :: h5dims(5)
      integer(HID_T) file_id

      FM=>FMap3D

      call h5open_f(h5err)
      call h5fopen_f('FootPrint.h5',H5F_ACC_RDONLY_F,file_id, h5err)

      h5dims(1)=1
      call h5ltread_dataset_f(file_id, "FMap3D/dim_a", H5T_NATIVE_INTEGER, FM%dim_a, h5dims(1:1), h5err)
      call h5ltread_dataset_f(file_id, "FMap3D/dim_gq", H5T_NATIVE_INTEGER, FM%dim_gq, h5dims(1:1), h5err)
      call h5ltread_dataset_f(file_id, "FMap3D/dim_gf", H5T_NATIVE_INTEGER, FM%dim_gf, h5dims(1:1), h5err)

      if(mype==0) write(stdout,'(3(A,I4))') "MG; FMap3D dim_a=",FMap3D%dim_a, "; dim_gq=",FMap3D%dim_gq, "; dim_gf=",FMap3D%dim_gf
      call flush(stdout)

      ! FM%dim_a  = 200
      ! FM%dim_gq = 300
      ! FM%dim_gf = 100

      Na = FM%dim_a
      Nq = FM%dim_gq
      Nf = FM%dim_gf

      allocate(FM%a_1D(Na),FM%gq_1D(Nq),FM%gf_1D(Nf))
      allocate(FM%Max_BB0_CCW(Nq,Na,Nf),FM%Max_BB0_CW(Nq,Na,Nf),FM%BB0(Nq,Na,Nf),FM%Mirror_Ratio(Nq,Na,Nf))

      h5dims(1)=Na
      call h5ltread_dataset_f(file_id, "FMap3D/a_1D", H5T_NATIVE_DOUBLE, FM%a_1D, h5dims(1:1), h5err)
      h5dims(1)=Nq
      call h5ltread_dataset_f(file_id, "FMap3D/gq_1D", H5T_NATIVE_DOUBLE, FM%gq_1D, h5dims(1:1), h5err)
      h5dims(1)=Nf
      call h5ltread_dataset_f(file_id, "FMap3D/gf_1D", H5T_NATIVE_DOUBLE, FM%gf_1D, h5dims(1:1), h5err)

      h5dims(1:3)=[Nq,Na,Nf]
      call h5ltread_dataset_f(file_id, "FMap3D/Max_BB0_CW", H5T_NATIVE_DOUBLE, FM%Max_BB0_CW, h5dims(1:3), h5err)
      call h5ltread_dataset_f(file_id, "FMap3D/Max_BB0_CCW", H5T_NATIVE_DOUBLE, FM%Max_BB0_CCW, h5dims(1:3), h5err)
      call h5ltread_dataset_f(file_id, "FMap3D/BB0", H5T_NATIVE_DOUBLE, FM%BB0, h5dims(1:3), h5err)
      call h5ltread_dataset_f(file_id, "FMap3D/M_eff", H5T_NATIVE_DOUBLE, FM%Mirror_Ratio, h5dims(1:3), h5err)
      call h5close_f(h5err)

      FM%inv_da = 1.0/(FM%a_1D(2)-FM%a_1D(1))
      FM%inv_dgq= 1.0/(FM%gq_1D(2)-FM%gq_1D(1))
      FM%inv_dgf= 1.0/(FM%gf_1D(2)-FM%gf_1D(1))


      if(mype==0) write(stdout,'(3(A,ES12.4))') "MG; FMap3D inv_da=",FMap3D%inv_da, "; inv_dgq=",FMap3D%inv_dgq, "; inv_dgf=",FMap3D%inv_dgf
      call flush(stdout)


   end subroutine Read_FMap3D_on_MG

   subroutine Read_ESI(filename)
      implicit none
      type(esi_data), pointer :: E
      character(128), intent(in) :: filename
      character(len=200) :: dummy
      integer :: Na,Nq,i,j
      real(wp) tmp(4)

      E=>esi

      open(unit=20201213,file=filename,action='read')

      read(20201213,*)
      read(20201213,*) E%dim_gq, dummy, E%dim_a
      read(20201213,*) dummy


      Na=E%dim_a
      Nq=E%dim_gq

      allocate(E%a_1D(Na),E%gq_1D(Nq))
      allocate(E%F(Na),E%Fa(Na),E%gFa(Na),E%gFaa(Na),E%gYa(Na),E%gYaa(Na),E%T(Na),E%Ta(Na),E%P(Na),E%Pa(Na))
      allocatE(E%R(Nq,Na),E%Ra(Nq,Na),E%Rq(Nq,Na),E%Raq(Nq,Na))
      allocate(E%Z(Nq,Na),E%Za(Nq,Na),E%Zq(Nq,Na),E%Zaq(Nq,Na))
      allocate(E%BB(Nq,Na),E%BBa(Nq,Na),E%BBq(Nq,Na),E%BBaq(Nq,Na))
      allocate(E%gh(Nq,Na),E%gh_a(Nq,Na),E%gh_q(Nq,Na),E%gh_aq(Nq,Na))
      allocate(E%ne(Na),E%dneda(Na))


      read(20201213,*) E%gq_1D(1:Nq)
      read(20201213,*)

      do i=1,Na
         read(20201213,*) E%a_1D(i), E%F(i), E%Fa(i)
      enddo
      read(20201213,*)

      do i=1,Na
         read(20201213,*) tmp(1:4)
         E%gFa(i)  = tmp(1)*E%a_1D(i)
         E%gFaa(i) = tmp(2)*E%a_1D(i) + tmp(1)
         E%gYa(i)  = tmp(3)*E%a_1D(i)
         E%gYaa(i) = tmp(4)*E%a_1D(i) + tmp(3)
      enddo
      read(20201213,*)


      do i=1,Na
         read(20201213,*) E%T(i), E%Ta(i), E%P(i), E%Pa(i)
      enddo
      read(20201213,*)


      do i=1,Na
         do j=1,Nq
            read(20201213,*) E%R(j,i), E%Ra(j,i), E%Rq(j,i), E%Raq(j,i)
         enddo
      enddo
      read(20201213,*)


      do i=1,Na
         do j=1,Nq
            read(20201213,*) E%Z(j,i), E%Za(j,i), E%Zq(j,i), E%Zaq(j,i)
         enddo
      enddo
      read(20201213,*)


      do i=1,Na
         do j=1,Nq
            read(20201213,*) E%BB(j,i), E%BBa(j,i), E%BBq(j,i), E%BBaq(j,i)
         enddo
      enddo
      read(20201213,*)

      do i=1,Na
         do j=1,Nq
            read(20201213,*) E%gh(j,i), E%gh_a(j,i), E%gh_q(j,i), E%gh_aq(j,i)
         enddo
      enddo
      read(20201213,*)

      do i=1,Na
         read(20201213,*) tmp(1:3)
         E%ne = 1e20*tmp(2)
         E%dneda = 1e20*tmp(3)
      enddo

      close(20201213)




      ! Post Processing


   end subroutine Read_ESI


   subroutine Write_ESI_hdf5
      use hdf5
      use h5lt
      implicit none

      character(len=128) filename
      integer(HID_T) file_id, gid
      integer :: h5err
      integer(HSIZE_T) :: h5dims(3)


      write(filename, '("ESI.h5")')

      call h5open_f(h5err)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)

      h5dims(1) = 1
      call h5ltmake_dataset_f(file_id, "dim_a",1,h5dims,H5T_NATIVE_INTEGER, esi%dim_a, h5err)
      call h5ltmake_dataset_f(file_id, "dim_gq",1,h5dims,H5T_NATIVE_INTEGER, esi%dim_gq, h5err)

      h5dims(1) = esi%dim_gq
      call h5ltmake_dataset_f(file_id, "gq_1D",1,h5dims,H5T_NATIVE_DOUBLE, esi%gq_1D, h5err)
      h5dims(1) = esi%dim_a
      call h5ltmake_dataset_f(file_id, "a_1D",1,h5dims,H5T_NATIVE_DOUBLE, esi%a_1D, h5err)
      call h5ltmake_dataset_f(file_id, "F",1,h5dims,H5T_NATIVE_DOUBLE, esi%F, h5err)
      call h5ltmake_dataset_f(file_id, "Fa",1,h5dims,H5T_NATIVE_DOUBLE, esi%Fa, h5err)
      call h5ltmake_dataset_f(file_id, "gFa",1,h5dims,H5T_NATIVE_DOUBLE, esi%gFa, h5err)
      call h5ltmake_dataset_f(file_id, "gFaa",1,h5dims,H5T_NATIVE_DOUBLE, esi%gFaa, h5err)
      call h5ltmake_dataset_f(file_id, "gYa",1,h5dims,H5T_NATIVE_DOUBLE, esi%gYa, h5err)
      call h5ltmake_dataset_f(file_id, "gYaa",1,h5dims,H5T_NATIVE_DOUBLE, esi%gYaa, h5err)
      call h5ltmake_dataset_f(file_id, "T",1,h5dims,H5T_NATIVE_DOUBLE, esi%T, h5err)
      call h5ltmake_dataset_f(file_id, "Ta",1,h5dims,H5T_NATIVE_DOUBLE, esi%Ta, h5err)
      call h5ltmake_dataset_f(file_id, "P",1,h5dims,H5T_NATIVE_DOUBLE, esi%P, h5err)
      call h5ltmake_dataset_f(file_id, "Pa",1,h5dims,H5T_NATIVE_DOUBLE, esi%Pa, h5err)


      ! 2D Data
      h5dims(1) = esi%dim_gq
      h5dims(2) = esi%dim_a
      call h5ltmake_dataset_f(file_id, "R",2,h5dims,H5T_NATIVE_DOUBLE, esi%R, h5err)
      call h5ltmake_dataset_f(file_id, "Ra",2,h5dims,H5T_NATIVE_DOUBLE, esi%Ra, h5err)
      call h5ltmake_dataset_f(file_id, "Rq",2,h5dims,H5T_NATIVE_DOUBLE, esi%Rq, h5err)
      call h5ltmake_dataset_f(file_id, "Raq",2,h5dims,H5T_NATIVE_DOUBLE, esi%Raq, h5err)
      call h5ltmake_dataset_f(file_id, "Z",2,h5dims,H5T_NATIVE_DOUBLE, esi%Z, h5err)
      call h5ltmake_dataset_f(file_id, "Za",2,h5dims,H5T_NATIVE_DOUBLE, esi%Za, h5err)
      call h5ltmake_dataset_f(file_id, "Zq",2,h5dims,H5T_NATIVE_DOUBLE, esi%Zq, h5err)
      call h5ltmake_dataset_f(file_id, "Zaq",2,h5dims,H5T_NATIVE_DOUBLE, esi%Zaq, h5err)
      call h5ltmake_dataset_f(file_id, "BB",2,h5dims,H5T_NATIVE_DOUBLE, esi%BB, h5err)
      call h5ltmake_dataset_f(file_id, "BBa",2,h5dims,H5T_NATIVE_DOUBLE, esi%BBa, h5err)
      call h5ltmake_dataset_f(file_id, "BBq",2,h5dims,H5T_NATIVE_DOUBLE, esi%BBq, h5err)
      call h5ltmake_dataset_f(file_id, "BBaq",2,h5dims,H5T_NATIVE_DOUBLE, esi%BBaq, h5err)
      call h5ltmake_dataset_f(file_id, "gh",2,h5dims,H5T_NATIVE_DOUBLE, esi%gh, h5err)
      call h5ltmake_dataset_f(file_id, "gh_a",2,h5dims,H5T_NATIVE_DOUBLE, esi%gh_a, h5err)
      call h5ltmake_dataset_f(file_id, "gh_q",2,h5dims,H5T_NATIVE_DOUBLE, esi%gh_q, h5err)
      call h5ltmake_dataset_f(file_id, "gh_aq",2,h5dims,H5T_NATIVE_DOUBLE, esi%gh_aq, h5err)
      call h5ltmake_dataset_f(file_id, "ne",2,h5dims,H5T_NATIVE_DOUBLE, esi%ne, h5err)
      call h5ltmake_dataset_f(file_id, "dneda",2,h5dims,H5T_NATIVE_DOUBLE, esi%dneda, h5err)

      call h5fclose_f(file_id, h5err)
      call h5close_f(h5err)


   end subroutine Write_ESI_hdf5

   subroutine Write_Magnetic_Grid_hdf5
      use hdf5
      use h5lt
      implicit none

      character(len=128) filename
      integer(HID_T) file_id, gid
      integer :: h5err
      integer(HSIZE_T) :: h5dims(5)


      write(filename, '("Magnetic_Grid.h5")')

      call h5open_f(h5err)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)

      h5dims(1) = 1
      call h5ltmake_dataset_f(file_id, "dim_a",1,h5dims,H5T_NATIVE_INTEGER, dim_a, h5err)
      call h5ltmake_dataset_f(file_id, "dim_gq",1,h5dims,H5T_NATIVE_INTEGER, dim_gq, h5err)
      call h5ltmake_dataset_f(file_id, "dim_gf",1,h5dims,H5T_NATIVE_INTEGER, dim_gf, h5err)


      h5dims(1) = dim_gf
      call h5ltmake_dataset_f(file_id, "gf_1D",1,h5dims,H5T_NATIVE_DOUBLE, gf_1D, h5err)
      h5dims(1) = dim_gq
      call h5ltmake_dataset_f(file_id, "gq_1D",1,h5dims,H5T_NATIVE_DOUBLE, gq_1D, h5err)
      h5dims(1) = dim_a
      call h5ltmake_dataset_f(file_id, "a_1D",1,h5dims,H5T_NATIVE_DOUBLE, a_1D, h5err)

      ! 2D Data
      h5dims(1) = dim_gq
      h5dims(2) = dim_a
      call h5ltmake_dataset_f(file_id, "R_2D",2,h5dims,H5T_NATIVE_DOUBLE, R_2D, h5err)
      call h5ltmake_dataset_f(file_id, "Z_2D",2,h5dims,H5T_NATIVE_DOUBLE, Z_2D, h5err)
      call h5ltmake_dataset_f(file_id, "Jacobian_2D",2,h5dims,H5T_NATIVE_DOUBLE, Jacobian_3D(:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "g_aa_cov_2D",2,h5dims,H5T_NATIVE_DOUBLE, g_aa_cov(:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "g_aq_cov_2D",2,h5dims,H5T_NATIVE_DOUBLE, g_aq_cov(:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "g_qq_cov_2D",2,h5dims,H5T_NATIVE_DOUBLE, g_qq_cov(:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "g_ff_cov_2D",2,h5dims,H5T_NATIVE_DOUBLE, g_ff_cov(:,:,1), h5err)


      h5dims(1:3)=[3,dim_gq,dim_a]
      call h5ltmake_dataset_f(file_id, "d12_R_2D",3,h5dims,H5T_NATIVE_DOUBLE, d12_R_2D, h5err)
      call h5ltmake_dataset_f(file_id, "d12_Z_2D",3,h5dims,H5T_NATIVE_DOUBLE, d12_Z_2D, h5err)
      call h5ltmake_dataset_f(file_id, "B0_cov",3,h5dims,H5T_NATIVE_DOUBLE, B0_cov(:,:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "B0_cntv",3,h5dims,H5T_NATIVE_DOUBLE, B0_cntv(:,:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "j0_cov",3,h5dims,H5T_NATIVE_DOUBLE, j0_cov(:,:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "j0_cntv",3,h5dims,H5T_NATIVE_DOUBLE, j0_cntv(:,:,:,1), h5err)


      h5dims(1:4)=[3,3,dim_gq,dim_a]
      call h5ltmake_dataset_f(file_id, "grad_B0_cov",4,h5dims,H5T_NATIVE_DOUBLE, grad_B0_cov(:,:,:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "grad_B0_cntv",4,h5dims,H5T_NATIVE_DOUBLE, grad_B0_cntv(:,:,:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "grad_j0_cov",4,h5dims,H5T_NATIVE_DOUBLE, grad_j0_cov(:,:,:,:,1), h5err)
      call h5ltmake_dataset_f(file_id, "grad_j0_cntv",4,h5dims,H5T_NATIVE_DOUBLE, grad_j0_cntv(:,:,:,:,1), h5err)


      h5dims(1:3)=[dim_gq,dim_a,dim_gf]
      call h5ltmake_dataset_f(file_id, "BB0",3,h5dims,H5T_NATIVE_DOUBLE, BB0, h5err)
      call h5ltmake_dataset_f(file_id, "dBB",3,h5dims,H5T_NATIVE_DOUBLE, dBB, h5err)
      call h5ltmake_dataset_f(file_id, "BB",3,h5dims,H5T_NATIVE_DOUBLE, BB, h5err)

      h5dims(1:4)=[3,dim_gq,dim_a,dim_gf]
      call h5ltmake_dataset_f(file_id, "grad_BB0",4,h5dims,H5T_NATIVE_DOUBLE, grad_BB0, h5err)
      call h5ltmake_dataset_f(file_id, "grad_dBB",4,h5dims,H5T_NATIVE_DOUBLE, grad_dBB, h5err)
      call h5ltmake_dataset_f(file_id, "grad_BB",4,h5dims,H5T_NATIVE_DOUBLE, grad_BB, h5err)


      h5dims(1:4)=[3,dim_gq,dim_a,dim_gf]
      call h5ltmake_dataset_f(file_id, "dB_cov",4,h5dims,H5T_NATIVE_DOUBLE, dB_cov, h5err)
      call h5ltmake_dataset_f(file_id, "dB_cntv",4,h5dims,H5T_NATIVE_DOUBLE, dB_cntv, h5err)
      call h5ltmake_dataset_f(file_id, "dj_cov",4,h5dims,H5T_NATIVE_DOUBLE, dj_cov, h5err)
      call h5ltmake_dataset_f(file_id, "dj_cntv",4,h5dims,H5T_NATIVE_DOUBLE, dj_cntv, h5err)

      h5dims(1:5)=[3,3,dim_gq,dim_a,dim_gf]
      call h5ltmake_dataset_f(file_id, "grad_dB_cov",5,h5dims,H5T_NATIVE_DOUBLE, grad_dB_cov, h5err)
      call h5ltmake_dataset_f(file_id, "grad_dB_cntv",5,h5dims,H5T_NATIVE_DOUBLE, grad_dB_cntv, h5err)
      call h5ltmake_dataset_f(file_id, "grad_dj_cov",5,h5dims,H5T_NATIVE_DOUBLE, grad_dj_cov, h5err)
      call h5ltmake_dataset_f(file_id, "grad_dj_cntv",5,h5dims,H5T_NATIVE_DOUBLE, grad_dj_cntv, h5err)

      h5dims(1:4)=[3,dim_gq,dim_a,dim_gf]
      call h5ltmake_dataset_f(file_id, "B_cov",4,h5dims,H5T_NATIVE_DOUBLE, B_cov, h5err)
      call h5ltmake_dataset_f(file_id, "B_cntv",4,h5dims,H5T_NATIVE_DOUBLE, B_cntv, h5err)
      call h5ltmake_dataset_f(file_id, "j_cov",4,h5dims,H5T_NATIVE_DOUBLE, j_cov, h5err)
      call h5ltmake_dataset_f(file_id, "j_cntv",4,h5dims,H5T_NATIVE_DOUBLE, j_cntv, h5err)

      h5dims(1:5)=[3,3,dim_gq,dim_a,dim_gf]
      call h5ltmake_dataset_f(file_id, "grad_B_cov",5,h5dims,H5T_NATIVE_DOUBLE, grad_B_cov, h5err)
      call h5ltmake_dataset_f(file_id, "grad_B_cntv",5,h5dims,H5T_NATIVE_DOUBLE, grad_B_cntv, h5err)
      call h5ltmake_dataset_f(file_id, "grad_j_cov",5,h5dims,H5T_NATIVE_DOUBLE, grad_j_cov, h5err)
      call h5ltmake_dataset_f(file_id, "grad_j_cntv",5,h5dims,H5T_NATIVE_DOUBLE, grad_j_cntv, h5err)

      h5dims(1:4)=[4,dim_gq,dim_a,dim_gf]
      call h5ltmake_dataset_f(file_id, "d23_BB0",4,h5dims,H5T_NATIVE_DOUBLE, d23_BB0, h5err)
      call h5ltmake_dataset_f(file_id, "d23_dBB",4,h5dims,H5T_NATIVE_DOUBLE, d23_dBB, h5err)
      call h5ltmake_dataset_f(file_id, "d23_BB",4,h5dims,H5T_NATIVE_DOUBLE, d23_BB, h5err)

      h5dims(1:5)=[4,3,dim_gq,dim_a,dim_gf]
      call h5ltmake_dataset_f(file_id, "d23_B_cov",5,h5dims,H5T_NATIVE_DOUBLE, d23_B_cov, h5err)


      call h5fclose_f(file_id, h5err)
      call h5close_f(h5err)


   end subroutine Write_Magnetic_Grid_hdf5


   subroutine Init_Magnetic_Grid
      use global_parameters, only:eq_data_name,pi,mype,stdout
      use mpi_defs
      implicit none

      integer i,j,k
      integer Na,Nq,Nf
      real(wp) a_min, a_max
      real(wp) tt0,tt1


      tt0=MPI_WTIME()

      call Read_ESI(eq_data_name)
      ! dim_a=SG%dim_a
      ! dim_gq=SG%dim_gq
      ! dim_gf=36
      dim_a=3*(esi%dim_a-1)+1   ! Finer than ESI but original nodes should be included
      dim_gq=3*(esi%dim_gq-1)+1
      dim_gf=65 ! 64 toroidal gaps that can deal with n=32 modes

      ! dim_a=esi%dim_a
      ! dim_gq=esi%dim_gq
      ! dim_gf=16


      Na=dim_a
      Nq=dim_gq
      Nf=dim_gf

      allocate(a_1D(dim_a),gq_1D(dim_gq),gf_1D(dim_gf))



      a_min=esi%a_1D(1)
      a_max=esi%a_1D(esi%dim_a)

      do i=1,dim_a
         a_1D(i) = a_min + real(i-1)*(a_max-a_min)/(dim_a-1)
      enddo

      do i=1,dim_gq
         gq_1D(i) = real(i-1)*2.0*pi/(dim_gq-1)
      enddo

      do i=1,dim_gf
         gf_1D(i) = real(i-1)*2.0*pi/(dim_gf-1)
      enddo

      da  = a_1D(2)-a_1D(1)
      dgq = gq_1D(2)-gq_1D(1)
      dgf = gf_1D(2)-gf_1D(1)

      inv_da  = 1.0/da
      inv_dgq = 1.0/dgq
      inv_dgf = 1.0/dgf


      allocate(R_2D(Nq,Na),Z_2D(Nq,Na))
      allocate(d12_R_2D(3,Nq,Na),d12_Z_2D(3,Nq,Na)) ! dRda, dRdgq, d2Rdagq
      allocate(Jacobian_3D(Nq,Na,Nf))
      allocate(g_aa_cov(Nq,Na,Nf),g_aq_cov(Nq,Na,Nf),g_qq_cov(Nq,Na,Nf),g_ff_cov(Nq,Na,Nf))

      ! Equilibrium B fields
      allocate(BB0(Nq,Na,Nf),grad_BB0(3,Nq,Na,Nf), d23_BB0(4,Nq,Na,Nf))
      allocate(B0_cov(3,Nq,Na,Nf), B0_cntv(3,Nq,Na,Nf), j0_cov(3,Nq,Na,Nf), j0_cntv(3,Nq,Na,Nf))
      allocate(grad_B0_cov(3,3,Nq,Na,Nf), grad_B0_cntv(3,3,Nq,Na,Nf), grad_j0_cov(3,3,Nq,Na,Nf), grad_j0_cntv(3,3,Nq,Na,Nf))
      allocate(d23_B0_cov(4,3,Nq,Na,Nf), d23_B0_cntv(4,3,Nq,Na,Nf), d23_j0_cov(4,3,Nq,Na,Nf), d23_j0_cntv(4,3,Nq,Na,Nf))

      ! dB (pertrubed) B fields
      allocate(dBB(Nq,Na,Nf),grad_dBB(3,Nq,Na,Nf), d23_dBB(4,Nq,Na,Nf))
      allocate(dB_cov(3,Nq,Na,Nf), dB_cntv(3,Nq,Na,Nf), dj_cov(3,Nq,Na,Nf), dj_cntv(3,Nq,Na,Nf))
      allocate(grad_dB_cov(3,3,Nq,Na,Nf), grad_dB_cntv(3,3,Nq,Na,Nf), grad_dj_cov(3,3,Nq,Na,Nf), grad_dj_cntv(3,3,Nq,Na,Nf))
      allocate(d23_dB_cov(4,3,Nq,Na,Nf), d23_dB_cntv(4,3,Nq,Na,Nf), d23_dj_cov(4,3,Nq,Na,Nf), d23_dj_cntv(4,3,Nq,Na,Nf))

      ! total B(=B0+dB)
      allocate(BB(Nq,Na,Nf),grad_BB(3,Nq,Na,Nf), d23_BB(4,Nq,Na,Nf))
      allocate(B_cov(3,Nq,Na,Nf),B_cntv(3,Nq,Na,Nf),j_cov(3,Nq,Na,Nf),j_cntv(3,Nq,Na,Nf))
      allocate(grad_B_cov(3,3,Nq,Na,Nf), grad_B_cntv(3,3,Nq,Na,Nf), grad_j_cov(3,3,Nq,Na,Nf), grad_j_cntv(3,3,Nq,Na,Nf))
      allocate(d23_B_cov(4,3,Nq,Na,Nf), d23_B_cntv(4,3,Nq,Na,Nf), d23_j_cov(4,3,Nq,Na,Nf), d23_j_cntv(4,3,Nq,Na,Nf))



      call Intp_from_ESI_to_MG
      call Cal_dB_from_Analytic_Alpha

      B_cov  = B0_cov  + dB_cov
      B_cntv = B0_cntv + dB_cntv
      j_cov  = j0_cov  + dj_cov
      j_cntv = j0_cntv + dj_cntv

      grad_B_cov  = grad_B0_cov  + grad_dB_cov
      grad_B_cntv = grad_B0_cntv + grad_dB_cntv
      grad_j_cov  = grad_j0_cov  + grad_dj_cov
      grad_j_cntv = grad_j0_cntv + grad_dj_cntv

      d23_B_cov  = d23_B0_cov  + d23_dB_cov
      d23_B_cntv = d23_B0_cntv + d23_dB_cntv
      d23_j_cov  = d23_j0_cov  + d23_dj_cov
      d23_j_cntv = d23_j0_cntv + d23_dj_cntv


      BB=sqrt(sum(B_cov*B_cntv,1))
      call Cal_grad_of_Scalar_on_MG(BB, grad_BB)
      call Cal_d23_of_Scalar_on_MG(BB,  d23_BB)

      if(mype==0) call Write_ESI_hdf5
      if(mype==0) call Write_Magnetic_Grid_hdf5
      ! if(mype==0) call test_MG_TriCubic_Interpolation

      !$acc update device(dim_a,dim_gq,dim_gf)
      !$acc update device(inv_da,inv_dgq,inv_dgf)
      !$acc update device(a_1D,gq_1D,gf_1D)
      !$acc update device(R_2D,d12_R_2D,Z_2D,d12_Z_2D)
      !$acc update device(BB,grad_BB,d23_BB)
      !$acc update device(B_cov,grad_B_cov,d23_B_cov)
      !$acc update device(B_cntv,grad_B_cntv,d23_B_cntv)
      !$acc update device(j_cntv,grad_j_cntv,d23_j_cntv)



      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "Init_Magnetic_Grid takes", tt1-tt0," secs"
      call flush(stdout)



   end subroutine Init_Magnetic_Grid

   subroutine test_MG_TriCubic_Interpolation
      use hdf5
      use h5lt
      use global_parameters
      use mpi_defs
      implicit none
      character(len=128) filename
      integer(HID_T) file_id, gid
      integer :: h5err
      integer(HSIZE_T) :: h5dims(4)
      real(wp),allocatable,dimension(:) :: test_a, test_gq, test_gf
      real(wp),allocatable,dimension(:,:) :: intp_R_2D
      real(wp),allocatable,dimension(:,:,:) :: intp_BB,intp_dBB, intp_grad_R_2D
      real(wp),allocatable,dimension(:,:,:,:) :: intp_grad_BB,intp_grad_dBB
      real(wp),allocatable,dimension(:,:,:,:) :: intp_B_cov, intp_B_cntv, intp_dB_cov, intp_dB_cntv
      integer test_Na, test_Nq, test_Nf
      integer i,j,k

      ! dim_gf=36
      ! dim_a=3*(esi%dim_a-1)+1   ! Finer than ESI but original nodes should be included
      ! dim_gq=3*(esi%dim_gq-1)+1
      ! dim_gf=50

      test_Na = 3*(dim_a-1)+1
      test_Nq = 3*(dim_gq-1)+1
      test_Nf = 3*(dim_gf-1)+1

      allocate(test_a(test_Na), test_gq(test_Nq), test_gf(test_Nf))
      allocate(intp_R_2D(test_Nq,test_Na), intp_grad_R_2D(2,test_Nq,test_Na))
      allocate(intp_BB(test_Nq,test_Na,test_Nf),intp_dBB(test_Nq,test_Na,test_Nf))
      allocate(intp_grad_BB(3,test_Nq,test_Na,test_Nf),intp_grad_dBB(3,test_Nq,test_Na,test_Nf))
      allocate(intp_B_cov(3,test_Nq,test_Na,test_Nf),intp_B_cntv(3,test_Nq,test_Na,test_Nf))
      allocate(intp_dB_cov(3,test_Nq,test_Na,test_Nf),intp_dB_cntv(3,test_Nq,test_Na,test_Nf))

      do i=1,test_Na
         test_a(i) = a_1D(1) + real(i-1)*(a_1D(dim_a)-a_1D(1))/(test_Na-1)
      enddo

      do i=1,test_Nq
         test_gq(i) = real(i-1)*2.0*pi/(test_Nq-1)
      enddo

      do i=1,test_Nf
         test_gf(i) = real(i-1)*2.0*pi/(test_Nf-1)
      enddo

      do i=1,test_Nq
         do j=1,test_Na
            call BiCubic_Interpolation_of_2D_Scalar_and_gradient(test_a(j),test_gq(i),R_2D, d12_R_2D,intp_R_2D(i,j), intp_grad_R_2D(:,i,j))
            do k=1,test_Nf
               ! call TriCubic_Interpolation_of_Scalar(test_a(j),test_gq(i),test_gf(k), BB, grad_BB, d23_BB, intp_BB(i,j,k))
               ! call TriCubic_Interpolation_of_Scalar(test_a(j),test_gq(i),test_gf(k), dBB, grad_dBB, d23_dBB, intp_dBB(i,j,k))
               call TriCubic_Interpolation_of_Scalar_and_gradient(test_a(j),test_gq(i),test_gf(k), BB, grad_BB, d23_BB, intp_BB(i,j,k), intp_grad_BB(:,i,j,k))
               call TriCubic_Interpolation_of_Scalar_and_gradient(test_a(j),test_gq(i),test_gf(k), dBB, grad_dBB, d23_dBB, intp_dBB(i,j,k), intp_grad_dBB(:,i,j,k))
               call TriCubic_Interpolation_of_Vec(test_a(j),test_gq(i),test_gf(k), B_cov, grad_B_cov, d23_B_cov, intp_B_cov(:,i,j,k))
               call TriCubic_Interpolation_of_Vec(test_a(j),test_gq(i),test_gf(k), B_cntv, grad_B_cntv, d23_B_cntv, intp_B_cntv(:,i,j,k))
               call TriCubic_Interpolation_of_Vec(test_a(j),test_gq(i),test_gf(k), dB_cov, grad_dB_cov, d23_dB_cov, intp_dB_cov(:,i,j,k))
               call TriCubic_Interpolation_of_Vec(test_a(j),test_gq(i),test_gf(k), dB_cntv, grad_dB_cntv, d23_dB_cntv, intp_dB_cntv(:,i,j,k))
            enddo
         enddo
      enddo


      write(filename, '("test_MG_TriCubicIntp.h5")')

      call h5open_f(h5err)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)
      h5dims(1)=test_Na
      call h5ltmake_dataset_f(file_id, "a_1D",1,h5dims,H5T_NATIVE_DOUBLE, test_a, h5err)
      h5dims(1)=test_Nq
      call h5ltmake_dataset_f(file_id, "gq_1D",1,h5dims,H5T_NATIVE_DOUBLE, test_gq, h5err)
      h5dims(1)=test_Nf
      call h5ltmake_dataset_f(file_id, "gf_1D",1,h5dims,H5T_NATIVE_DOUBLE, test_gf, h5err)

      h5dims(1:2) = [test_Nq,test_Na]
      call h5ltmake_dataset_f(file_id, "R_2D",2,h5dims,H5T_NATIVE_DOUBLE, intp_R_2D, h5err)
      h5dims(1:3) = [2,test_Nq,test_Na]
      call h5ltmake_dataset_f(file_id, "grad_R_2D",3,h5dims,H5T_NATIVE_DOUBLE, intp_grad_R_2D, h5err)

      h5dims(1:3) = [test_Nq,test_Na,test_Nf]
      call h5ltmake_dataset_f(file_id, "BB",3,h5dims,H5T_NATIVE_DOUBLE, intp_BB, h5err)
      call h5ltmake_dataset_f(file_id, "dBB",3,h5dims,H5T_NATIVE_DOUBLE, intp_dBB, h5err)
      h5dims(1:4) = [3,test_Nq,test_Na,test_Nf]
      call h5ltmake_dataset_f(file_id, "grad_BB",4,h5dims,H5T_NATIVE_DOUBLE, intp_grad_BB, h5err)
      call h5ltmake_dataset_f(file_id, "grad_dBB",4,h5dims,H5T_NATIVE_DOUBLE, intp_grad_dBB, h5err)
      call h5ltmake_dataset_f(file_id, "B_cov",4,h5dims,H5T_NATIVE_DOUBLE, intp_B_cov, h5err)
      call h5ltmake_dataset_f(file_id, "B_cntv",4,h5dims,H5T_NATIVE_DOUBLE, intp_B_cntv, h5err)
      call h5ltmake_dataset_f(file_id, "dB_cntv",4,h5dims,H5T_NATIVE_DOUBLE, intp_dB_cntv, h5err)
      call h5ltmake_dataset_f(file_id, "dB_cov",4,h5dims,H5T_NATIVE_DOUBLE, intp_dB_cov, h5err)

      call h5fclose_f(file_id, h5err)
      call h5close_f(h5err)


      deallocate(test_a,test_gq,test_gf)
      deallocate(intp_BB,intp_dBB,intp_B_cov,intp_B_cntv)
      deallocate(intp_dB_cntv,intp_dB_cov)
      deallocate(intp_R_2D, intp_grad_R_2D)
      deallocate(intp_grad_BB,intp_grad_dBB)

   end subroutine test_MG_TriCubic_Interpolation





   subroutine Intp_from_ESI_to_MG
      implicit none

      real(wp) ha00,ha10,ha01,ha11, dha00,dha10,dha01,dha11
      real(wp) hq00,hq10,hq01,hq11, dhq00,dhq10,dhq01,dhq11
      integer aid,gqid
      real(wp) a,gq
      real(wp) esi_da,inv_esi_da,esi_dgq,inv_esi_dgq
      real(wp) at,qt
      integer esi_aid,esi_gqid
      real(wp) Jacob, esi_N, esi_K
      real(wp) F,dFda,gYa,P,T
      real(wp) R,dRda,dRdq,d2Rdadq,Z,dZda,dZdq,d2Zdadq
      real(wp) f1D(0:1), dfa1D(0:1), f2D(0:1,0:1), dfa2D(0:1,0:1), dfq2D(0:1,0:1), dfaq2D(0:1,0:1)
      real(wp) f_Ha(0:1), dfa_Ha(0:1), dfq_Ha(0:1), dfaq_Ha(0:1) ! values at Half a


      esi_da =esi%a_1D(2)-esi%a_1D(1)
      esi_dgq=esi%gq_1D(2)-esi%gq_1D(1)

      inv_esi_da  = 1.0/(esi_da)
      inv_esi_dgq = 1.0/(esi_dgq)

      do aid=1,dim_a
         a = a_1D(aid)
         do gqid=1,dim_gq
            gq=gq_1D(gqid)

            esi_aid  = min(esi%dim_a-1,floor((a-esi%a_1D(1))*inv_esi_da) +1)
            esi_gqid = floor((gq-esi%gq_1D(1))*inv_esi_dgq) +1

            at=(a-esi%a_1D(esi_aid))*inv_esi_da
            qt=(gq-esi%gq_1D(esi_gqid))*inv_esi_dgq

            ha00 = 2*at**3 -3*at**2 +1
            ha10 = (at**3 -2*at**2 +at)*esi_da
            ha01 = -2*at**3 +3*at**2
            ha11 = (at**3 -at**2)*esi_da
            dha00 = (6*at**2 -6*at)*inv_esi_da
            dha10 = 3*at**2 -4*at +1
            dha01 = (-6*at**2 +6*at)*inv_esi_da
            dha11 = 3*at**2 -2*at

            hq00 = 2*qt**3 -3*qt**2 +1
            hq10 = (qt**3 -2*qt**2 +qt)*esi_dgq
            hq01 = -2*qt**3 +3*qt**2
            hq11 = (qt**3 -qt**2)*esi_dgq
            dhq00 = (6*qt**2 -6*qt)*inv_esi_dgq
            dhq10 = 3*qt**2 -4*qt +1
            dhq01 = (-6*qt**2 +6*qt)*inv_esi_dgq
            dhq11 = 3*qt**2 -2*qt


            ! Interpolate 1D quantities first
            f1D  = esi%F(esi_aid:esi_aid+1)  ! function value
            dfa1D= esi%Fa(esi_aid:esi_aid+1) ! derivative
            F    =  ha00*f1D(0) + ha10*dfa1D(0)  + ha01*f1D(1) +  ha11*dfa1D(1)
            dFda = dha00*f1D(0) +dha10*dfa1D(0)  +dha01*f1D(1) + dha11*dfa1D(1)

            f1D  = esi%gYa(esi_aid:esi_aid+1)  ! function value
            dfa1D= esi%gYaa(esi_aid:esi_aid+1) ! derivative
            gYa  = ha00*f1D(0) + ha10*dfa1D(0)  + ha01*f1D(1) +  ha11*dfa1D(1)

            f1D  = esi%P(esi_aid:esi_aid+1)  ! function value
            dfa1D= esi%Pa(esi_aid:esi_aid+1) ! derivative
            P    = ha00*f1D(0) + ha10*dfa1D(0)  + ha01*f1D(1) +  ha11*dfa1D(1)

            f1D  = esi%T(esi_aid:esi_aid+1)  ! function value
            dfa1D= esi%Ta(esi_aid:esi_aid+1) ! derivative
            T    = ha00*f1D(0) + ha10*dfa1D(0)  + ha01*f1D(1) +  ha11*dfa1D(1)

            ! Interpolate 2D quantatities
            f2D    = esi%R(esi_gqid:esi_gqid+1,esi_aid:esi_aid+1)
            dfa2D  = esi%Ra(esi_gqid:esi_gqid+1,esi_aid:esi_aid+1)
            dfq2D  = esi%Rq(esi_gqid:esi_gqid+1,esi_aid:esi_aid+1)
            dfaq2D = esi%Raq(esi_gqid:esi_gqid+1,esi_aid:esi_aid+1)
            f_Ha   =  ha00*f2D(:,0)   + ha10*dfa2D(:,0)  + ha01*f2D(:,1)   + ha11*dfa2D(:,1) ! f value at a Half (q,a)=(0,0.5)&(1,0.5)
            dfa_Ha = dha00*f2D(:,0)   +dha10*dfa2D(:,0)  +dha01*f2D(:,1)   +dha11*dfa2D(:,1)
            dfq_Ha =  ha00*dfq2D(:,0) + ha10*dfaq2D(:,0) + ha01*dfq2D(:,1) + ha11*dfaq2D(:,1)
            dfaq_Ha= dha00*dfq2D(:,0) +dha10*dfaq2D(:,0) +dha01*dfq2D(:,1) +dha11*dfaq2D(:,1)
            R      =  hq00*f_Ha(0)    + hq10*dfq_Ha(0)   + hq01*f_Ha(1)    + hq11*dfq_Ha(1)
            dRdq   = dhq00*f_Ha(0)    +dhq10*dfq_Ha(0)   +dhq01*f_Ha(1)    +dhq11*dfq_Ha(1)
            dRda   =  hq00*dfa_Ha(0)  + hq10*dfaq_Ha(0)  + hq01*dfa_Ha(1)  + hq11*dfaq_Ha(1)
            d2Rdadq= dhq00*dfa_Ha(0)  +dhq10*dfaq_Ha(0)  +dhq01*dfa_Ha(1)  +dhq11*dfaq_Ha(1)

            f2D    = esi%Z(esi_gqid:esi_gqid+1,esi_aid:esi_aid+1)
            dfa2D  = esi%Za(esi_gqid:esi_gqid+1,esi_aid:esi_aid+1)
            dfq2D  = esi%Zq(esi_gqid:esi_gqid+1,esi_aid:esi_aid+1)
            dfaq2D = esi%Zaq(esi_gqid:esi_gqid+1,esi_aid:esi_aid+1)
            f_Ha   =  ha00*f2D(:,0)   + ha10*dfa2D(:,0)  + ha01*f2D(:,1)   + ha11*dfa2D(:,1) ! f value at a Half (q,a)=(0,0.5)&(1,0.5)
            dfa_Ha = dha00*f2D(:,0)   +dha10*dfa2D(:,0)  +dha01*f2D(:,1)   +dha11*dfa2D(:,1)
            dfq_Ha =  ha00*dfq2D(:,0) + ha10*dfaq2D(:,0) + ha01*dfq2D(:,1) + ha11*dfaq2D(:,1)
            dfaq_Ha= dha00*dfq2D(:,0) +dha10*dfaq2D(:,0) +dha01*dfq2D(:,1) +dha11*dfaq2D(:,1)
            Z      =  hq00*f_Ha(0)    + hq10*dfq_Ha(0)   + hq01*f_Ha(1)    + hq11*dfq_Ha(1)
            dZdq   = dhq00*f_Ha(0)    +dhq10*dfq_Ha(0)   +dhq01*f_Ha(1)    +dhq11*dfq_Ha(1)
            dZda   =  hq00*dfa_Ha(0)  + hq10*dfaq_Ha(0)  + hq01*dfa_Ha(1)  + hq11*dfaq_Ha(1)
            d2Zdadq= dhq00*dfa_Ha(0)  +dhq10*dfaq_Ha(0)  +dhq01*dfa_Ha(1)  +dhq11*dfaq_Ha(1)

            Jacob  = R*(dZda*dRdq-dZdq*dRda)


            R_2D(gqid,aid) = R
            Z_2D(gqid,aid) = Z
            d12_R_2D(:,gqid,aid) = [dRda, dRdq, d2Rdadq]
            d12_Z_2D(:,gqid,aid) = [dZda, dZdq, d2Zdadq]

            Jacobian_3D(gqid,aid,:) = Jacob


            esi_N = (dRda*dRdq+dZda*dZdq)/Jacob
            esi_K = (dRdq**2+dZdq**2)/Jacob

            ! B0_cov(:,gqid,aid,:)  = [-esi_N*gYa, -esi_K*gYa, F]
            ! B0_cntv(:,gqid,aid,:) = [0.0, -gYa/Jacob, F/(R*R)]
            ! j0_cntv(:,gqid,aid,:) = [0.0, -dFda/Jacob, P+T/(R*R)]
            B0_cov(1,gqid,aid,:)  = -esi_N*gYa
            B0_cov(2,gqid,aid,:)  = -esi_K*gYa
            B0_cov(3,gqid,aid,:)  = F

            B0_cntv(1,gqid,aid,:) = 0.0
            B0_cntv(2,gqid,aid,:) = -gYa/Jacob
            B0_cntv(3,gqid,aid,:) = F/(R*R)

            j0_cntv(1,gqid,aid,:) = 0.0
            j0_cntv(2,gqid,aid,:) = -dFda/Jacob
            j0_cntv(3,gqid,aid,:) = P+T/(R*R)


            g_aa_cov(gqid,aid,:) = dRda**2+dZda**2
            g_aq_cov(gqid,aid,:) = dRda*dRdq+dZda*dZdq
            g_qq_cov(gqid,aid,:) = dRdq**2+dZdq**2
            g_ff_cov(gqid,aid,:) = R*R

            j0_cov(1,gqid,aid,:) = g_aq_cov(gqid,aid,:)*j0_cntv(2,gqid,aid,:)
            j0_cov(2,gqid,aid,:) = g_qq_cov(gqid,aid,:)*j0_cntv(2,gqid,aid,:)
            j0_cov(3,gqid,aid,:) = g_ff_cov(gqid,aid,:)*j0_cntv(3,gqid,aid,:)
         enddo
      enddo

      BB0=sqrt(sum(B0_cov*B0_cntv,1))
      call Cal_grad_of_Scalar_on_MG(BB0, grad_BB0)
      call Cal_d23_of_Scalar_on_MG(BB0,  d23_BB0)

      call Cal_grad_of_Vec_on_MG(B0_cov, grad_B0_cov)
      call Cal_grad_of_Vec_on_MG(B0_cntv,grad_B0_cntv)
      call Cal_grad_of_Vec_on_MG(j0_cov, grad_j0_cov)
      call Cal_grad_of_Vec_on_MG(j0_cntv,grad_j0_cntv)

      call Cal_d23_of_Vec_on_MG(B0_cov,  d23_B0_cov)
      call Cal_d23_of_Vec_on_MG(B0_cntv, d23_B0_cntv)
      call Cal_d23_of_Vec_on_MG(j0_cov,  d23_j0_cov)
      call Cal_d23_of_Vec_on_MG(j0_cntv, d23_j0_cntv)



   end subroutine Intp_from_ESI_to_MG


   subroutine Cal_dB_from_Analytic_Alpha
      implicit none

      real(wp) a,gq,gf
      real(wp) gaE,grad_ga_cov(3),gatE
      integer aa,qq,ff


      do ff=1,dim_gf
         gf=gf_1D(ff)
         do aa=1,dim_a
            a=a_1D(aa)
            do qq=1,dim_gq
               gq=gq_1D(qq)

               call Get_MP_alpha_amplitude(a,gq,gf,0.0,gaE,grad_ga_cov,gatE)
               dB_cntv(1,qq,aa,ff) = gaE*j0_cntv(1,qq,aa,ff) + (grad_ga_cov(2)*B0_cov(3,qq,aa,ff)-grad_ga_cov(3)*B0_cov(2,qq,aa,ff))/Jacobian_3D(qq,aa,ff)
               dB_cntv(2,qq,aa,ff) = gaE*j0_cntv(2,qq,aa,ff) + (grad_ga_cov(3)*B0_cov(1,qq,aa,ff)-grad_ga_cov(1)*B0_cov(3,qq,aa,ff))/Jacobian_3D(qq,aa,ff)
               dB_cntv(3,qq,aa,ff) = gaE*j0_cntv(3,qq,aa,ff) + (grad_ga_cov(1)*B0_cov(2,qq,aa,ff)-grad_ga_cov(2)*B0_cov(1,qq,aa,ff))/Jacobian_3D(qq,aa,ff)
            enddo
         enddo
      enddo


      dB_cov(1,:,:,:) = g_aa_cov*dB_cntv(1,:,:,:) + g_aq_cov*dB_cntv(2,:,:,:)
      dB_cov(2,:,:,:) = g_aq_cov*dB_cntv(1,:,:,:) + g_qq_cov*dB_cntv(2,:,:,:)
      dB_cov(3,:,:,:) = g_ff_cov*dB_cntv(3,:,:,:)

      call Cal_grad_of_Vec_on_MG(dB_cntv,grad_dB_cntv)
      call Cal_grad_of_Vec_on_MG(dB_cov, grad_dB_cov)

      dj_cntv(1,:,:,:) = (grad_dB_cov(2,3,:,:,:)-grad_dB_cov(3,2,:,:,:))/Jacobian_3D
      dj_cntv(2,:,:,:) = (grad_dB_cov(3,1,:,:,:)-grad_dB_cov(1,3,:,:,:))/Jacobian_3D
      dj_cntv(3,:,:,:) = (grad_dB_cov(1,2,:,:,:)-grad_dB_cov(2,1,:,:,:))/Jacobian_3D

      dj_cov(1,:,:,:) = g_aa_cov*dj_cntv(1,:,:,:) + g_aq_cov*dj_cntv(2,:,:,:)
      dj_cov(2,:,:,:) = g_aq_cov*dj_cntv(1,:,:,:) + g_qq_cov*dj_cntv(2,:,:,:)
      dj_cov(3,:,:,:) = g_ff_cov*dj_cntv(3,:,:,:)


      call Cal_grad_of_Vec_on_MG(dj_cntv,grad_dj_cntv)
      call Cal_grad_of_Vec_on_MG(dj_cov, grad_dj_cov)

      call Cal_d23_of_Vec_on_MG(dB_cov,  d23_dB_cov)
      call Cal_d23_of_Vec_on_MG(dB_cntv, d23_dB_cntv)
      call Cal_d23_of_Vec_on_MG(dj_cov,  d23_dj_cov)
      call Cal_d23_of_Vec_on_MG(dj_cntv, d23_dj_cntv)

      dBB=sqrt(sum(dB_cov*dB_cntv,1))
      call Cal_grad_of_Scalar_on_MG(dBB, grad_dBB)
      call Cal_d23_of_Scalar_on_MG(dBB,  d23_dBB)




   end subroutine Cal_dB_from_Analytic_Alpha


   subroutine Cal_grad_of_Scalar_on_MG(Sca,gradSca)
      implicit none
      real(wp), intent(in) :: Sca(dim_gq,dim_a,dim_gf)
      real(wp), intent(out):: gradSca(3,dim_gq,dim_a,dim_gf)
      integer Nq,Na,Nf

      Na=dim_a
      Nq=dim_gq
      Nf=dim_gf

      ! dSca/da
      gradSca(1,:,1,:) = inv_da*(Sca(:,2,:)-Sca(:,1,:))
      gradSca(1,:,2:Na-1,:) = 0.5*inv_da*(Sca(:,3:Na,:)-Sca(:,1:Na-2,:))
      gradSca(1,:,Na,:) =inv_da*(Sca(:,Na,:)-Sca(:,Na-1,:))

      ! dSca/dgq
      gradSca(2,1,:,:) = 0.5*inv_dgq*(Sca(2,:,:)-Sca(Nq-1,:,:))
      gradSca(2,2:Nq-1,:,:) = 0.5*inv_dgq*(Sca(3:Nq,:,:)-Sca(1:Nq-2,:,:))
      gradSca(2,Nq,:,:) = gradSca(2,1,:,:)

      ! dSca/dgf
      gradSca(3,:,:,1) = 0.5*inv_dgf*(Sca(:,:,2)-Sca(:,:,Nf-1))
      gradSca(3,:,:,2:Nf-1) = 0.5*inv_dgf*(Sca(:,:,3:Nf)-Sca(:,:,1:Nf-2))
      gradSca(3,:,:,Nf) = gradSca(3,:,:,1)

   end subroutine Cal_grad_of_Scalar_on_MG

   subroutine Cal_grad_of_Vec_on_MG(Vec,gradVec)
      implicit none
      real(wp), intent(in) :: Vec(3,dim_gq,dim_a,dim_gf)
      real(wp), intent(out):: gradVec(3,3,dim_gq,dim_a,dim_gf)
      integer Nq,Na,Nf

      Na=dim_a
      Nq=dim_gq
      Nf=dim_gf

      ! dVec/da
      gradVec(1,:,:,1,:) = inv_da*(Vec(:,:,2,:)-Vec(:,:,1,:))
      gradVec(1,:,:,2:Na-1,:) = 0.5*inv_da*(Vec(:,:,3:Na,:)-Vec(:,:,1:Na-2,:))
      gradVec(1,:,:,Na,:) =inv_da*(Vec(:,:,Na,:)-Vec(:,:,Na-1,:))

      ! dVec/dgq
      gradVec(2,:,1,:,:) = 0.5*inv_dgq*(Vec(:,2,:,:)-Vec(:,Nq-1,:,:))
      gradVec(2,:,2:Nq-1,:,:) = 0.5*inv_dgq*(Vec(:,3:Nq,:,:)-Vec(:,1:Nq-2,:,:))
      gradVec(2,:,Nq,:,:) = gradVec(2,:,1,:,:)

      ! dVec/dgf
      gradVec(3,:,:,:,1) = 0.5*inv_dgf*(Vec(:,:,:,2)-Vec(:,:,:,Nf-1))
      gradVec(3,:,:,:,2:Nf-1) = 0.5*inv_dgf*(Vec(:,:,:,3:Nf)-Vec(:,:,:,1:Nf-2))
      gradVec(3,:,:,:,Nf) = gradVec(3,:,:,:,1)


   end subroutine Cal_grad_of_Vec_on_MG


   subroutine Cal_d23_of_Scalar_on_MG(in_Sca,out_d23Sca)
      implicit none
      real(wp), intent(in) :: in_Sca(dim_gq,dim_a,dim_gf)
      real(wp), intent(out):: out_d23Sca(4,dim_gq,dim_a,dim_gf)
      real(wp) Sca(0:dim_gq+1,dim_a,0:dim_gf+1)
      integer Nq,Na,Nf

      Na=dim_a
      Nq=dim_gq
      Nf=dim_gf

      ! Copy in_Sca into Sca, and build ghost cells for periodicity in (gq,gf) direction
      Sca(1:Nq,:,1:Nf) = in_Sca
      Sca(0,:,:) = Sca(Nq-1,:,:)
      Sca(Nq+1,:,:) = Sca(2,:,:)
      Sca(:,:,0) = Sca(:,:,Nf-1)
      Sca(:,:,Nf+1) = Sca(:,:,2)

      out_d23Sca=0.0

      !d23Sca(1,:,:,:)  = d^2(Sca)/(da*dgq)
      out_d23Sca(1,:,2:Na-1,:) = 0.25*inv_da*inv_dgq*(Sca(2:Nq+1,3:Na,1:Nf)+Sca(0:Nq-1,1:Na-2,1:Nf)-Sca(2:Nq+1,1:Na-2,1:Nf)-Sca(0:Nq-1,3:Na,1:Nf))

      !d23Sca(2,:,:,:)  = d^2(Sca)/(da*dgf)
      out_d23Sca(2,:,2:Na-1,:) = 0.25*inv_da*inv_dgf*(Sca(1:Nq,3:Na,2:Nf+1)+Sca(1:Nq,1:Na-2,0:Nf-1)-Sca(1:Nq,3:Na,0:Nf-1)-Sca(1:Nq,1:Na-2,2:Nf+1))

      !d23Sca(3,:,:,:)  = d^2(Sca)/(dgq*dgf)
      out_d23Sca(3,:,:,:) = 0.25*inv_dgq*inv_dgf*(Sca(2:Nq+1,1:Na,2:Nf+1)+Sca(0:Nq-1,1:Na,0:Nf-1)-Sca(2:Nq+1,1:Na,0:Nf-1)-Sca(0:Nq-1,1:Na,2:Nf+1))

      !d23Sca(4,:,:,:)  = d^3(Sca)/(da*dgq*dgf)
      out_d23Sca(4,:,2:Na-1,:) = 0.125*inv_da*inv_dgq*inv_dgf &
                                       *(Sca(2:Nq+1,3:Na,2:Nf+1)+Sca(0:Nq-1,1:Na-2,2:Nf+1)+Sca(0:Nq-1,3:Na,0:Nf-1)+Sca(2:Nq+1,1:Na-2,0:Nf-1) &
                                        -Sca(0:Nq-1,3:Na,2:Nf+1)-Sca(2:Nq+1,1:Na-2,2:Nf+1)-Sca(2:Nq+1,3:Na,0:Nf-1)-Sca(0:Nq-1,1:Na-2,0:Nf-1))

   end subroutine Cal_d23_of_Scalar_on_MG


   subroutine Cal_d23_of_Vec_on_MG(in_Vec,out_d23Vec)
      implicit none
      real(wp), intent(in) :: in_Vec(3,dim_gq,dim_a,dim_gf)
      real(wp), intent(out):: out_d23Vec(4,3,dim_gq,dim_a,dim_gf)
      real(wp) Vec(3,0:dim_gq+1,dim_a,0:dim_gf+1)
      integer Nq,Na,Nf

      Na=dim_a
      Nq=dim_gq
      Nf=dim_gf

      ! Copy in_Vec into Vec, and build ghost cells for periodicity in (gq,gf) direction
      Vec(:,1:Nq,:,1:Nf) = in_Vec
      Vec(:,0,:,:) = Vec(:,Nq-1,:,:)
      Vec(:,Nq+1,:,:) = Vec(:,2,:,:)
      Vec(:,:,:,0) = Vec(:,:,:,Nf-1)
      Vec(:,:,:,Nf+1) = Vec(:,:,:,2)

      out_d23Vec=0.0

      !d23Vec(1,:,:,:,:)  = d^2(Vec)/(da*dgq)
      out_d23Vec(1,:,:,2:Na-1,:) = 0.25*inv_da*inv_dgq*(Vec(:,2:Nq+1,3:Na,1:Nf)+Vec(:,0:Nq-1,1:Na-2,1:Nf)-Vec(:,2:Nq+1,1:Na-2,1:Nf)-Vec(:,0:Nq-1,3:Na,1:Nf))

      !d23Vec(2,:,:,:,:)  = d^2(Vec)/(da*dgf)
      out_d23Vec(2,:,:,2:Na-1,:) = 0.25*inv_da*inv_dgf*(Vec(:,1:Nq,3:Na,2:Nf+1)+Vec(:,1:Nq,1:Na-2,0:Nf-1)-Vec(:,1:Nq,3:Na,0:Nf-1)-Vec(:,1:Nq,1:Na-2,2:Nf+1))

      !d23Vec(3,:,:,:,:)  = d^2(Vec)/(dgq*dgf)
      out_d23Vec(3,:,:,:,:) = 0.25*inv_dgq*inv_dgf*(Vec(:,2:Nq+1,1:Na,2:Nf+1)+Vec(:,0:Nq-1,1:Na,0:Nf-1)-Vec(:,2:Nq+1,1:Na,0:Nf-1)-Vec(:,0:Nq-1,1:Na,2:Nf+1))

      !d23Vec(4,:,:,:,:)  = d^3(Vec)/(da*dgq*dgf)
      out_d23Vec(4,:,:,2:Na-1,:) = 0.125*inv_da*inv_dgq*inv_dgf &
                                       *(Vec(:,2:Nq+1,3:Na,2:Nf+1)+Vec(:,0:Nq-1,1:Na-2,2:Nf+1)+Vec(:,0:Nq-1,3:Na,0:Nf-1)+Vec(:,2:Nq+1,1:Na-2,0:Nf-1) &
                                        -Vec(:,0:Nq-1,3:Na,2:Nf+1)-Vec(:,2:Nq+1,1:Na-2,2:Nf+1)-Vec(:,2:Nq+1,3:Na,0:Nf-1)-Vec(:,0:Nq-1,1:Na-2,0:Nf-1))



   end subroutine Cal_d23_of_Vec_on_MG


   ! subroutine Cal_d23_of_Scalar_on_MG(Sca,d23Sca)
   !    implicit none
   !    real(wp), intent(in) :: Sca(dim_gq,dim_a,dim_gf)
   !    real(wp), intent(out):: d23Sca(4,dim_gq,dim_a,dim_gf)
   !    integer Nq,Na,Nf

   !    Na=dim_a
   !    Nq=dim_gq
   !    Nf=dim_gf

   !    !d23Sca(1,:,:,:)  = d^2(Sca)/(da*dgq)
   !    d23Sca(1,1,2:Na-1,:) = 0.25*inv_da*inv_dgq*(Sca(2,3:Na,:)+Sca(Nq-1,1:Na-2,:)-Sca(2,1:Na-2,:)-Sca(Nq-1,3:Na,:))
   !    d23Sca(1,2:Nq-1,2:Na-1,:) = 0.25*inv_da*inv_dgq*(Sca(3:Nq,3:Na,:)+Sca(1:Nq-2,1:Na-2,:)-Sca(3:Nq,1:Na-2,:)-Sca(1:Nq-2,3:Na,:))
   !    d23Sca(1,Nq,2:Na-1,:) = d23Sca(1,1,2:Na-1,:)


   !    !d23Sca(2,:,:,:)  = d^2(Sca)/(da*dgf)
   !    d23Sca(2,:,2:Na-1,1) = 0.25*inv_da*inv_dgf*(Sca(:,3:Na,2)+Sca(:,1:Na-2,Nf-1)-Sca(:,3:Na,Nf-1)-Sca(:,1:Na-2,2))
   !    d23Sca(2,:,2:Na-1,2:Nf-1) = 0.25*inv_da*inv_dgf*(Sca(:,3:Na,3:Nf)+Sca(:,1:Na-2,1:Nf-2)-Sca(:,3:Na,1:Nf-2)-Sca(:,1:Na-2,3:Nf))
   !    d23Sca(2,:,2:Na-1,Nf) = d23Sca(2,:,2:Na-1,1)


   !    !d23Sca(3,:,:,:)  = d^2(Sca)/(dgq*dgf)
   !    d23Sca(3,1,:,2:Nf-1) = 0.25*inv_dgq*inv_dgf*(Sca(2,:,3:Nf)+Sca(Nq-1,:,1:Nf-2)-Sca(2,:,1:Nf-2)-Sca(Nq-1,:,3:Nf))
   !    d23Sca(3,2:Nq-1,:,1) = 0.25*inv_dgq*inv_dgf*(Sca(3:Nq,:,2)+Sca(1:Nq-2,:,Nf-1)-Sca(3:Nq,:,Nf-1)-Sca(1:Nq-2,:,2))
   !    d23Sca(3,2:Nq-1,:,2:Nf-1) = 0.25*inv_dgq*inv_dgf*(Sca(3:Nq,:,3:Nf)+Sca(1:Nq-2,:,1:Nf-2)-Sca(3:Nq,:,1:Nf-2)-Sca(1:Nq-2,:,3:Nf))
   !    d23Sca(3,Nq,:,2:Nf-1) = d23Sca(3,1,:,2:Nf-1)
   !    d23Sca(3,2:Nq-1,:,Nf) = d23Sca(3,2:Nq-1,:,1)


   !    !d23Sca(4,:,:,:)  = d^3(Sca)/(da*dgq*dgf)
   !    d23Sca(4,1,2:Na-1,2:Nf-1) = 0.125*inv_da*inv_dgq*inv_dgf &
   !                                     *(Sca(2,3:Na,3:Nf)+Sca(Nq-1,1:Na-2,3:Nf)+Sca(Nq-1,3:Na,1:Nf-2)+Sca(2,1:Na-2,1:Nf-2) &
   !                                      -Sca(Nq-1,3:Na,3:Nf)-Sca(2,1:Na-2,3:Nf)-Sca(2,3:Na,1:Nf-2)-Sca(Nq-1,1:Na-2,1:Nf-2))
   !    d23Sca(4,2:Nq-1,2:Na-1,1) = 0.125*inv_da*inv_dgq*inv_dgf &
   !                                     *(Sca(3:Nq,3:Na,2)+Sca(1:Nq-2,1:Na-2,2)+Sca(1:Nq-2,3:Na,Nf-1)+Sca(3:Nq,1:Na-2,Nf-1) &
   !                                      -Sca(1:Nq-2,3:Na,2)-Sca(3:Nq,1:Na-2,2)-Sca(3:Nq,3:Na,Nf-1)-Sca(1:Nq-2,1:Na-2,Nf-1))
   !    d23Sca(4,2:Nq-1,2:Na-1,2:Nf-1) = 0.125*inv_da*inv_dgq*inv_dgf &
   !                                     *(Sca(3:Nq,3:Na,3:Nf)+Sca(1:Nq-2,1:Na-2,3:Nf)+Sca(1:Nq-2,3:Na,1:Nf-2)+Sca(3:Nq,1:Na-2,1:Nf-2) &
   !                                      -Sca(1:Nq-2,3:Na,3:Nf)-Sca(3:Nq,1:Na-2,3:Nf)-Sca(3:Nq,3:Na,1:Nf-2)-Sca(1:Nq-2,1:Na-2,1:Nf-2))
   !    d23Sca(4,Nq,2:Na-1,2:Nf-1) = d23Sca(4,1,2:Na-1,2:Nf-1)
   !    d23Sca(4,2:Nq-1,2:Na-1,Nf) = d23Sca(4,2:Nq-1,2:Na-1,1)

   ! end subroutine Cal_d23_of_Scalar_on_MG



   ! subroutine Cal_d23_of_Vec_on_MG(Vec,d23Vec)
   !    implicit none
   !    real(wp), intent(in) :: Vec(3,dim_gq,dim_a,dim_gf)
   !    real(wp), intent(out):: d23Vec(4,3,dim_gq,dim_a,dim_gf)
   !    integer Nq,Na,Nf

   !    Na=dim_a
   !    Nq=dim_gq
   !    Nf=dim_gf

   !    ! !d23Vec(1,:,:,:,:)  = d^2(Vec)/(da*dgq)
   !    ! !d23Vec(2,:,:,:,:)  = d^2(Vec)/(da*dgf)
   !    ! !d23Vec(3,:,:,:,:)  = d^2(Vec)/(dgq*dgf)
   !    ! d23Vec(1,:,2:Nq-1,2:Na-1,:) = 0.25*inv_da*inv_dgq*(Vec(:,3:Nq,3:Na,:)+Vec(:,1:Nq-2,1:Na-2,:)-Vec(:,3:Nq,1:Na-2,:)-Vec(:,1:Nq-2,3:Na,:))
   !    ! d23Vec(2,:,:,2:Na-1,2:Nf-1) = 0.25*inv_da*inv_dgf*(Vec(:,:,3:Na,3:Nf)+Vec(:,:,1:Na-2,1:Nf-2)-Vec(:,:,3:Na,1:Nf-2)-Vec(:,:,1:Na-2,3:Nf))
   !    ! d23Vec(3,:,2:Nq-1,:,2:Nf-1) = 0.25*inv_dgq*inv_dgf*(Vec(:,3:Nq,:,3:Nf)+Vec(:,1:Nq-2,:,1:Nf-2)-Vec(:,3:Nq,:,1:Nf-2)-Vec(:,1:Nq-2,:,3:Nf))

   !    ! !d23Vec(4,:,:,:,:)  = d^3(Vec)/(da*dgq*dgf)
   !    ! d23Vec(4,:,2:Nq-1,2:Na-1,2:Nf-1) = 0.125*inv_da*inv_dgq*inv_dgf &
   !    !                                     *(Vec(:,3:Nq,3:Na,3:Nf)+Vec(:,1:Nq-2,1:Na-2,3:Nf)+Vec(:,1:Nq-2,3:Na,1:Nf-2)+Vec(:,3:Nq,1:Na-2,1:Nf-2) &
   !    !                                      -Vec(:,1:Nq-2,3:Na,3:Nf)-Vec(:,3:Nq,1:Na-2,3:Nf)-Vec(:,3:Nq,3:Na,1:Nf-2)-Vec(:,1:Nq-2,1:Na-2,1:Nf-2))


   !    !d23Vec(1,:,:,:,:)  = d^2(Vec)/(da*dgq)
   !    d23Vec(1,:,1,2:Na-1,:) = 0.25*inv_da*inv_dgq*(Vec(:,2,3:Na,:)+Vec(:,Nq-1,1:Na-2,:)-Vec(:,2,1:Na-2,:)-Vec(:,Nq-1,3:Na,:))
   !    d23Vec(1,:,2:Nq-1,2:Na-1,:) = 0.25*inv_da*inv_dgq*(Vec(:,3:Nq,3:Na,:)+Vec(:,1:Nq-2,1:Na-2,:)-Vec(:,3:Nq,1:Na-2,:)-Vec(:,1:Nq-2,3:Na,:))
   !    d23Vec(1,:,Nq,2:Na-1,:) = d23Vec(1,:,1,2:Na-1,:)


   !    !d23Vec(2,:,:,:,:)  = d^2(Vec)/(da*dgf)
   !    d23Vec(2,:,:,2:Na-1,1) = 0.25*inv_da*inv_dgf*(Vec(:,:,3:Na,2)+Vec(:,:,1:Na-2,Nf-1)-Vec(:,:,3:Na,Nf-1)-Vec(:,:,1:Na-2,2))
   !    d23Vec(2,:,:,2:Na-1,2:Nf-1) = 0.25*inv_da*inv_dgf*(Vec(:,:,3:Na,3:Nf)+Vec(:,:,1:Na-2,1:Nf-2)-Vec(:,:,3:Na,1:Nf-2)-Vec(:,:,1:Na-2,3:Nf))
   !    d23Vec(2,:,:,2:Na-1,Nf) = d23Vec(2,:,:,2:Na-1,1)


   !    !d23Vec(3,:,:,:,:)  = d^2(Vec)/(dgq*dgf)
   !    d23Vec(3,:,1,:,2:Nf-1) = 0.25*inv_dgq*inv_dgf*(Vec(:,2,:,3:Nf)+Vec(:,Nq-1,:,1:Nf-2)-Vec(:,2,:,1:Nf-2)-Vec(:,Nq-1,:,3:Nf))
   !    d23Vec(3,:,2:Nq-1,:,1) = 0.25*inv_dgq*inv_dgf*(Vec(:,3:Nq,:,2)+Vec(:,1:Nq-2,:,Nf-1)-Vec(:,3:Nq,:,Nf-1)-Vec(:,1:Nq-2,:,2))
   !    d23Vec(3,:,2:Nq-1,:,2:Nf-1) = 0.25*inv_dgq*inv_dgf*(Vec(:,3:Nq,:,3:Nf)+Vec(:,1:Nq-2,:,1:Nf-2)-Vec(:,3:Nq,:,1:Nf-2)-Vec(:,1:Nq-2,:,3:Nf))
   !    d23Vec(3,:,Nq,:,2:Nf-1) = d23Vec(3,:,1,:,2:Nf-1)
   !    d23Vec(3,:,2:Nq-1,:,Nf) = d23Vec(3,:,2:Nq-1,:,1)


   !    !d23Vec(4,:,:,:,:)  = d^3(Vec)/(da*dgq*dgf)
   !    d23Vec(4,:,1,2:Na-1,2:Nf-1) = 0.125*inv_da*inv_dgq*inv_dgf &
   !                                     *(Vec(:,2,3:Na,3:Nf)+Vec(:,Nq-1,1:Na-2,3:Nf)+Vec(:,Nq-1,3:Na,1:Nf-2)+Vec(:,2,1:Na-2,1:Nf-2) &
   !                                      -Vec(:,Nq-1,3:Na,3:Nf)-Vec(:,2,1:Na-2,3:Nf)-Vec(:,2,3:Na,1:Nf-2)-Vec(:,Nq-1,1:Na-2,1:Nf-2))
   !    d23Vec(4,:,2:Nq-1,2:Na-1,1) = 0.125*inv_da*inv_dgq*inv_dgf &
   !                                     *(Vec(:,3:Nq,3:Na,2)+Vec(:,1:Nq-2,1:Na-2,2)+Vec(:,1:Nq-2,3:Na,Nf-1)+Vec(:,3:Nq,1:Na-2,Nf-1) &
   !                                      -Vec(:,1:Nq-2,3:Na,2)-Vec(:,3:Nq,1:Na-2,2)-Vec(:,3:Nq,3:Na,Nf-1)-Vec(:,1:Nq-2,1:Na-2,Nf-1))
   !    d23Vec(4,:,2:Nq-1,2:Na-1,2:Nf-1) = 0.125*inv_da*inv_dgq*inv_dgf &
   !                                     *(Vec(:,3:Nq,3:Na,3:Nf)+Vec(:,1:Nq-2,1:Na-2,3:Nf)+Vec(:,1:Nq-2,3:Na,1:Nf-2)+Vec(:,3:Nq,1:Na-2,1:Nf-2) &
   !                                      -Vec(:,1:Nq-2,3:Na,3:Nf)-Vec(:,3:Nq,1:Na-2,3:Nf)-Vec(:,3:Nq,3:Na,1:Nf-2)-Vec(:,1:Nq-2,1:Na-2,1:Nf-2))
   !    d23Vec(4,:,Nq,2:Na-1,2:Nf-1) = d23Vec(4,:,1,2:Na-1,2:Nf-1)
   !    d23Vec(4,:,2:Nq-1,2:Na-1,Nf) = d23Vec(4,:,2:Nq-1,2:Na-1,1)


   ! end subroutine Cal_d23_of_Vec_on_MG


   subroutine BiCubic_Interpolation_of_2D_Scalar(in_a, in_gq, in_Sca, in_d12Sca, out_Sca)
      !$acc routine seq
      implicit none
      real(wp),intent(in) :: in_a,in_gq
      real(wp),intent(in) :: in_Sca(dim_gq,dim_a)
      real(wp),intent(in) :: in_d12Sca(3,dim_gq,dim_a) ! First dimension is for three elements of dFda, dFdgq, d2Fdadgq
      real(wp),intent(out) :: out_Sca
      real(wp) h00(2),h10(2),h01(2),h11(2)
      integer aid,qid
      real(wp) at,qt,tt(2),du(2)
      real(wp) F(0:1,0:1),d12F(3,0:1,0:1) ! Values at corners (2D in a,q)
      real(wp) F_Ha(0:1),dFq_Ha(0:1) ! Values at Half a (1D in q)


      aid = max(1,min(dim_a-1, floor((in_a-a_1D(1))*inv_da +1)))
      qid = max(1,min(dim_gq-1, floor((in_gq-gq_1D(1))*inv_dgq +1)))

      at = (in_a-a_1D(aid))*inv_da
      qt = (in_gq-gq_1D(qid))*inv_dgq

      tt = [at,qt]
      du = [da,dgq]

      h00 = 2*tt**3 -3*tt**2 +1
      h10 = (tt**3 -2*tt**2 +tt)*du
      h01 = -2*tt**3 +3*tt**2
      h11 = (tt**3 -tt**2)*du

      F  = in_Sca(qid:qid+1,aid:aid+1)
      d12F = in_d12Sca(:,qid:qid+1,aid:aid+1)

      ! Find values at Half a (for 1D in gq)
      F_Ha = h00(1)*F(:,0) +h10(1)*d12F(1,:,0) +h01(1)*F(:,1) +h11(1)*d12F(1,:,1)
      dFq_Ha = h00(1)*d12F(2,:,0) +h10(1)*d12F(3,:,0) +h01(1)*d12F(2,:,1) +h11(1)*d12F(3,:,1)

      ! Find values at Half a & Half gq (at the given point)
      out_Sca = h00(2)*F_Ha(0) +h10(2)*dFq_Ha(0) +h01(2)*F_Ha(1) +h11(2)*dFq_Ha(1) ! F

   end subroutine BiCubic_Interpolation_of_2D_Scalar

   subroutine BiCubic_Interpolation_of_2D_Scalar_and_gradient(in_a, in_gq, in_Sca, in_d12Sca, out_Sca, out_gradSca)
      !$acc routine seq
      implicit none
      real(wp),intent(in) :: in_a,in_gq
      real(wp),intent(in) :: in_Sca(dim_gq,dim_a)
      real(wp),intent(in) :: in_d12Sca(3,dim_gq,dim_a) ! First dimension is for three elements of dFda, dFdgq, d2Fdadgq
      real(wp),intent(out) :: out_Sca, out_gradSca(2)
      real(wp) h00(2),h10(2),h01(2),h11(2)
      real(wp) dh00(2),dh10(2),dh01(2),dh11(2)
      integer aid,qid
      real(wp) at,qt,tt(2),du(2)
      real(wp) F(0:1,0:1),d12F(3,0:1,0:1) ! Values at corners (2D in a,q)
      real(wp) F_Ha(0:1),dFa_Ha(0:1),dFq_Ha(0:1),dFaq_Ha(0:1) ! Values at Half a (1D in q)


      aid = max(1,min(dim_a-1, floor((in_a-a_1D(1))*inv_da +1)))
      qid = max(1,min(dim_gq-1, floor((in_gq-gq_1D(1))*inv_dgq +1)))

      at = (in_a-a_1D(aid))*inv_da
      qt = (in_gq-gq_1D(qid))*inv_dgq


      tt = [at,qt]
      du = [da,dgq]

      h00 = 2*tt**3 -3*tt**2 +1
      h10 = (tt**3 -2*tt**2 +tt)*du
      h01 = -2*tt**3 +3*tt**2
      h11 = (tt**3 -tt**2)*du
      dh00 = (6*tt**2 -6*tt)/du
      dh10 = 3*tt**2 -4*tt +1
      dh01 = (-6*tt**2 +6*tt)/du
      dh11 = 3*tt**2 -2*tt

      F  = in_Sca(qid:qid+1,aid:aid+1)
      d12F = in_d12Sca(:,qid:qid+1,aid:aid+1)

      ! Find values at Half a (for 1D in gq)
      F_Ha = h00(1)*F(:,0) +h10(1)*d12F(1,:,0) +h01(1)*F(:,1) +h11(1)*d12F(1,:,1)
      dFa_Ha = dh00(1)*F(:,0) +dh10(1)*d12F(1,:,0) +dh01(1)*F(:,1) +dh11(1)*d12F(1,:,1)
      dFq_Ha = h00(1)*d12F(2,:,0) +h10(1)*d12F(3,:,0) +h01(1)*d12F(2,:,1) +h11(1)*d12F(3,:,1)
      dFaq_Ha = dh00(1)*d12F(2,:,0) +dh10(1)*d12F(3,:,0) +dh01(1)*d12F(2,:,1) +dh11(1)*d12F(3,:,1)

      ! Find values at Half a & Half gq (at the given point)
      out_Sca = h00(2)*F_Ha(0) +h10(2)*dFq_Ha(0) +h01(2)*F_Ha(1) +h11(2)*dFq_Ha(1) ! F
      out_gradSca(1) = h00(2)*dFa_Ha(0) +h10(2)*dFaq_Ha(0) +h01(2)*dFa_Ha(1) +h11(2)*dFaq_Ha(1) ! dFda
      out_gradSca(2) = dh00(2)*F_Ha(0) +dh10(2)*dFq_Ha(0) +dh01(2)*F_Ha(1) +dh11(2)*dFq_Ha(1) ! dFdgq

   end subroutine BiCubic_Interpolation_of_2D_Scalar_and_gradient


   subroutine TriCubic_Interpolation_of_Scalar(in_a,in_gq,in_gf, in_Sca, in_gradSca, in_d23Sca, out_Sca)
      !$acc routine seq
      implicit none
      real(wp),intent(in) :: in_a,in_gq,in_gf
      real(wp),intent(in) :: in_Sca(dim_gq,dim_a,dim_gf)
      real(wp),intent(in) :: in_gradSca(3,dim_gq,dim_a,dim_gf), in_d23Sca(4,dim_gq,dim_a,dim_gf)
      real(wp),intent(out) :: out_Sca
      real(wp) h00(3),h10(3),h01(3),h11(3)
      ! real(wp) dh00(3),dh10(3),dh01(3),dh11(3)
      integer aid,qid,fid
      real(wp) at,qt,ft,tt(3),du(3)
      real(wp) F(0:1,0:1,0:1),dF(3,0:1,0:1,0:1),d23F(4,0:1,0:1,0:1) ! Values at corners (3D in a,q,f)
      real(wp) F_Ha(0:1,0:1), dFq_Ha(0:1,0:1), dFf_Ha(0:1,0:1), d2Fqf_Ha(0:1,0:1)  ! Values at Half a (2D in q,f)
      real(wp) F_Haq(0:1), dFf_Haq(0:1) ! Values at Half a & Half q (1D in f)


      aid = max(1,min(dim_a-1, floor((in_a-a_1D(1))*inv_da +1)))
      qid = max(1,min(dim_gq-1, floor((in_gq-gq_1D(1))*inv_dgq +1)))
      fid = max(1,min(dim_gf-1, floor((in_gf-gf_1D(1))*inv_dgf +1)))


      at = (in_a-a_1D(aid))*inv_da
      qt = (in_gq-gq_1D(qid))*inv_dgq
      ft = (in_gf-gf_1D(fid))*inv_dgf


      tt = [at,qt,ft]
      du = [da,dgq,dgf]

      h00 = 2*tt**3 -3*tt**2 +1
      h10 = (tt**3 -2*tt**2 +tt)*du
      h01 = -2*tt**3 +3*tt**2
      h11 = (tt**3 -tt**2)*du
      ! dh00 = (6*tt**2 -6*tt)/du
      ! dh10 = 3*tt**2 -4*tt +1
      ! dh01 = (-6*tt**2 +6*tt)/du
      ! dh11 = 3*tt**2 -2*tt

      F  = in_Sca(qid:qid+1,aid:aid+1,fid:fid+1)
      dF = in_gradSca(:,qid:qid+1,aid:aid+1,fid:fid+1)
      d23F = in_d23Sca(:,qid:qid+1,aid:aid+1,fid:fid+1)


      ! Find values at Half a (for 2D in gq,gf)
      F_Ha = h00(1)*F(:,0,:) +h10(1)*dF(1,:,0,:) +h01(1)*F(:,1,:) +h11(1)*dF(1,:,1,:)
      dFq_Ha = h00(1)*dF(2,:,0,:) +h10(1)*d23F(1,:,0,:) +h01(1)*dF(2,:,1,:) +h11(1)*d23F(1,:,1,:)
      dFf_Ha = h00(1)*dF(3,:,0,:) +h10(1)*d23F(2,:,0,:) +h01(1)*dF(3,:,1,:) +h11(1)*d23F(2,:,1,:)
      d2Fqf_Ha = h00(1)*d23F(3,:,0,:) +h10(1)*d23F(4,:,0,:) +h01(1)*d23F(3,:,1,:) +h11(1)*d23F(4,:,1,:)

      ! Find values at Half a & Half gq (for 1D in gf)
      F_Haq = h00(2)*F_Ha(0,:) +h10(2)*dFq_Ha(0,:) +h01(2)*F_Ha(1,:) +h11(2)*dFq_Ha(1,:)
      dFf_Haq =  h00(2)*dFf_Ha(0,:) +h10(2)*d2Fqf_Ha(0,:) +h01(2)*dFf_Ha(1,:) +h11(2)*d2Fqf_Ha(1,:)

      ! Find interpolated value at the given point
      out_Sca = h00(3)*F_Haq(0) +h10(3)*dFf_Haq(0) +h01(3)*F_Haq(1) +h11(3)*dFf_Haq(1)

   end subroutine TriCubic_Interpolation_of_Scalar

   subroutine TriCubic_Interpolation_of_Scalar_and_gradient(in_a,in_gq,in_gf, in_Sca, in_gradSca, in_d23Sca, out_Sca, out_gradSca)
      !$acc routine seq
      implicit none
      real(wp),intent(in) :: in_a,in_gq,in_gf
      real(wp),intent(in) :: in_Sca(dim_gq,dim_a,dim_gf)
      real(wp),intent(in) :: in_gradSca(3,dim_gq,dim_a,dim_gf), in_d23Sca(4,dim_gq,dim_a,dim_gf)
      real(wp),intent(out) :: out_Sca
      real(wp),intent(out) :: out_gradSca(3)
      real(wp) h00(3),h10(3),h01(3),h11(3)
      real(wp) dh00(3),dh10(3),dh01(3),dh11(3)
      integer aid,qid,fid
      real(wp) at,qt,ft,tt(3),du(3)
      real(wp) F(0:1,0:1,0:1),dF(3,0:1,0:1,0:1),d23F(4,0:1,0:1,0:1) ! Values at corners (3D in a,q,f)
      real(wp) F_Ha(0:1,0:1), dFa_Ha(0:1,0:1), dFq_Ha(0:1,0:1), dFf_Ha(0:1,0:1), d2Faq_Ha(0:1,0:1), d2Faf_Ha(0:1,0:1), d2Fqf_Ha(0:1,0:1), d3Faqf_Ha(0:1,0:1)  ! Values at Half a (2D in q,f)
      real(wp) F_Haq(0:1), dFa_Haq(0:1), dFq_Haq(0:1), dFf_Haq(0:1), d2Faf_Haq(0:1), d2Fqf_Haq(0:1) ! Values at Half a & Half q (1D in f)


      aid = max(1,min(dim_a-1, floor((in_a-a_1D(1))*inv_da +1)))
      qid = max(1,min(dim_gq-1, floor((in_gq-gq_1D(1))*inv_dgq +1)))
      fid = max(1,min(dim_gf-1, floor((in_gf-gf_1D(1))*inv_dgf +1)))


      at = (in_a-a_1D(aid))*inv_da
      qt = (in_gq-gq_1D(qid))*inv_dgq
      ft = (in_gf-gf_1D(fid))*inv_dgf


      tt = [at,qt,ft]
      du = [da,dgq,dgf]

      h00 = 2*tt**3 -3*tt**2 +1
      h10 = (tt**3 -2*tt**2 +tt)*du
      h01 = -2*tt**3 +3*tt**2
      h11 = (tt**3 -tt**2)*du
      dh00 = (6*tt**2 -6*tt)/du
      dh10 = 3*tt**2 -4*tt +1
      dh01 = (-6*tt**2 +6*tt)/du
      dh11 = 3*tt**2 -2*tt

      F  = in_Sca(qid:qid+1,aid:aid+1,fid:fid+1)
      dF = in_gradSca(:,qid:qid+1,aid:aid+1,fid:fid+1)
      d23F = in_d23Sca(:,qid:qid+1,aid:aid+1,fid:fid+1)


      ! Find values at Half a (for 2D in gq,gf)
      F_Ha = h00(1)*F(:,0,:) +h10(1)*dF(1,:,0,:) +h01(1)*F(:,1,:) +h11(1)*dF(1,:,1,:)
      dFa_Ha=dh00(1)*F(:,0,:) +dh10(1)*dF(1,:,0,:) +dh01(1)*F(:,1,:) +dh11(1)*dF(1,:,1,:)
      dFq_Ha = h00(1)*dF(2,:,0,:) +h10(1)*d23F(1,:,0,:) +h01(1)*dF(2,:,1,:) +h11(1)*d23F(1,:,1,:)
      dFf_Ha = h00(1)*dF(3,:,0,:) +h10(1)*d23F(2,:,0,:) +h01(1)*dF(3,:,1,:) +h11(1)*d23F(2,:,1,:)
      d2Faq_Ha = dh00(1)*dF(2,:,0,:) +dh10(1)*d23F(1,:,0,:) +dh01(1)*dF(2,:,1,:) +dh11(1)*d23F(1,:,1,:)
      d2Faf_Ha = dh00(1)*dF(3,:,0,:) +dh10(1)*d23F(2,:,0,:) +dh01(1)*dF(3,:,1,:) +dh11(1)*d23F(2,:,1,:)
      d2Fqf_Ha = h00(1)*d23F(3,:,0,:) +h10(1)*d23F(4,:,0,:) +h01(1)*d23F(3,:,1,:) +h11(1)*d23F(4,:,1,:)
      d3Faqf_Ha=dh00(1)*d23F(3,:,0,:) +dh10(1)*d23F(4,:,0,:) +dh01(1)*d23F(3,:,1,:) +dh11(1)*d23F(4,:,1,:)

      ! Find values at Half a & Half gq (for 1D in gf)
      F_Haq = h00(2)*F_Ha(0,:) +h10(2)*dFq_Ha(0,:) +h01(2)*F_Ha(1,:) +h11(2)*dFq_Ha(1,:)
      dFa_Haq =  h00(2)*dFa_Ha(0,:) +h10(2)*d2Faq_Ha(0,:) +h01(2)*dFa_Ha(0,:) +h11(2)*d2Faq_Ha(0,:)
      dFq_Haq = dh00(2)*F_Ha(0,:) +dh10(2)*dFq_Ha(0,:) +dh01(2)*F_Ha(1,:) +dh11(2)*dFq_Ha(1,:)
      dFf_Haq =  h00(2)*dFf_Ha(0,:) +h10(2)*d2Fqf_Ha(0,:) +h01(2)*dFf_Ha(1,:) +h11(2)*d2Fqf_Ha(1,:)
      d2Faf_Haq = h00(2)*d2Faf_Ha(0,:) +h10(2)*d3Faqf_Ha(0,:) +h01(2)*d2Faf_Ha(1,:) +h11(2)*d3Faqf_Ha(1,:)
      d2Fqf_Haq = dh00(2)*dFf_Ha(0,:) +dh10(2)*d2Fqf_Ha(0,:) +dh01(2)*dFf_Ha(1,:) +dh11(2)*d2Fqf_Ha(1,:)



      ! Find interpolated value at the given point
      out_Sca = h00(3)*F_Haq(0) +h10(3)*dFf_Haq(0) +h01(3)*F_Haq(1) +h11(3)*dFf_Haq(1)
      out_gradSca(1) = h00(3)*dFa_Haq(0) +h10(3)*d2Faf_Haq(0) +h01(3)*dFa_Haq(1) +h11(3)*d2Faf_Haq(1)
      out_gradSca(2) = h00(3)*dFq_Haq(0) +h10(3)*d2Fqf_Haq(0) +h01(3)*dFq_Haq(1) +h11(3)*d2Fqf_Haq(1)
      out_gradSca(3) = dh00(3)*F_Haq(0) +dh10(3)*dFf_Haq(0) +dh01(3)*F_Haq(1) +dh11(3)*dFf_Haq(1)

   end subroutine TriCubic_Interpolation_of_Scalar_and_gradient


   subroutine TriCubic_Interpolation_of_Vec(in_a,in_gq,in_gf, in_Vec, in_gradVec, in_d23Vec, out_Vec)
      !$acc routine seq
      implicit none
      real(wp),intent(in) :: in_a,in_gq,in_gf
      real(wp),intent(in) :: in_Vec(3,dim_gq,dim_a,dim_gf)
      real(wp),intent(in) :: in_gradVec(3,3,dim_gq,dim_a,dim_gf), in_d23Vec(4,3,dim_gq,dim_a,dim_gf)
      real(wp),intent(out) :: out_Vec(3)
      real(wp) h00(3),h10(3),h01(3),h11(3)
      ! real(wp) dh00(3),dh10(3),dh01(3),dh11(3)
      integer aid,qid,fid
      real(wp) at,qt,ft,tt(3),du(3)
      real(wp) F(3,0:1,0:1,0:1),dF(3,3,0:1,0:1,0:1),d23F(4,3,0:1,0:1,0:1) ! Values at corners (3D in a,q,f)
      real(wp) F_Ha(3,0:1,0:1), dFq_Ha(3,0:1,0:1), dFf_Ha(3,0:1,0:1), d2Fqf_Ha(3,0:1,0:1)  ! Values at Half a (2D in q,f)
      real(wp) F_Haq(3,0:1), dFf_Haq(3,0:1) ! Values at Half a & Half q (1D in f)


      aid = max(1,min(dim_a-1, floor((in_a-a_1D(1))*inv_da +1)))
      qid = max(1,min(dim_gq-1, floor((in_gq-gq_1D(1))*inv_dgq +1)))
      fid = max(1,min(dim_gf-1, floor((in_gf-gf_1D(1))*inv_dgf +1)))

      at = (in_a-a_1D(aid))*inv_da
      qt = (in_gq-gq_1D(qid))*inv_dgq
      ft = (in_gf-gf_1D(fid))*inv_dgf


      tt = [at,qt,ft]
      du = [da,dgq,dgf]

      h00 = 2*tt**3 -3*tt**2 +1
      h10 = (tt**3 -2*tt**2 +tt)*du
      h01 = -2*tt**3 +3*tt**2
      h11 = (tt**3 -tt**2)*du
      ! dh00 = (6*tt**2 -6*tt)/du
      ! dh10 = 3*tt**2 -4*tt +1
      ! dh01 = (-6*tt**2 +6*tt)/du
      ! dh11 = 3*tt**2 -2*tt


      F  = in_Vec(:,qid:qid+1,aid:aid+1,fid:fid+1)
      dF = in_gradVec(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
      d23F = in_d23Vec(:,:,qid:qid+1,aid:aid+1,fid:fid+1)


      ! Find values at Half a (for 2D in gq,gf)
      F_Ha = h00(1)*F(:,:,0,:) +h10(1)*dF(1,:,:,0,:) +h01(1)*F(:,:,1,:) +h11(1)*dF(1,:,:,1,:)
      dFq_Ha = h00(1)*dF(2,:,:,0,:) +h10(1)*d23F(1,:,:,0,:) +h01(1)*dF(2,:,:,1,:) +h11(1)*d23F(1,:,:,1,:)
      dFf_Ha = h00(1)*dF(3,:,:,0,:) +h10(1)*d23F(2,:,:,0,:) +h01(1)*dF(3,:,:,1,:) +h11(1)*d23F(2,:,:,1,:)
      d2Fqf_Ha = h00(1)*d23F(3,:,:,0,:) +h10(1)*d23F(4,:,:,0,:) +h01(1)*d23F(3,:,:,1,:) +h11(1)*d23F(4,:,:,1,:)

      ! Find values at Half a & Half gq (for 1D in gf)
      F_Haq = h00(2)*F_Ha(:,0,:) +h10(2)*dFq_Ha(:,0,:) +h01(2)*F_Ha(:,1,:) +h11(2)*dFq_Ha(:,1,:)
      dFf_Haq =  h00(2)*dFf_Ha(:,0,:) +h10(2)*d2Fqf_Ha(:,0,:) +h01(2)*dFf_Ha(:,1,:) +h11(2)*d2Fqf_Ha(:,1,:)

      ! Find interpolated value at the given point
      out_Vec = h00(3)*F_Haq(:,0) +h10(3)*dFf_Haq(:,0) +h01(3)*F_Haq(:,1) +h11(3)*dFf_Haq(:,1)

   end subroutine TriCubic_Interpolation_of_Vec


   subroutine Cubic_Interpolation_of_VariousStuff_for_GCmotion(in_a,in_gq,in_gf,P_R,P_grad_R,P_Z,P_grad_Z,P_BB,P_grad_BB,P_B_cov,P_B_cntv,P_j_cntv)
      !$acc routine seq
      implicit none
      ! real(wp),intent(in) :: in_a,in_gq,in_gf
      ! real(wp),intent(in) :: in_Vec(3,dim_gq,dim_a,dim_gf)
      ! real(wp),intent(in) :: in_gradVec(3,3,dim_gq,dim_a,dim_gf), in_d23Vec(4,3,dim_gq,dim_a,dim_gf)
      ! real(wp),intent(out) :: out_Vec(3)
      real(wp) h00(3),h10(3),h01(3),h11(3)
      real(wp) dh00(3),dh10(3),dh01(3),dh11(3)
      integer aid,qid,fid
      real(wp) at,qt,ft,tt(3),du(3)
      real(wp) S(0:1,0:1,0:1),dS(3,0:1,0:1,0:1),d23S(4,0:1,0:1,0:1) ! Values at corners (3D in a,q,f)
      real(wp) S_Ha(0:1,0:1), dSa_Ha(0:1,0:1), dSq_Ha(0:1,0:1), dSf_Ha(0:1,0:1), d2Saq_Ha(0:1,0:1), d2Saf_Ha(0:1,0:1), d2Sqf_Ha(0:1,0:1), d3Saqf_Ha(0:1,0:1)  ! Values at Half a (2D in q,f)
      real(wp) S_Haq(0:1), dSa_Haq(0:1), dSq_Haq(0:1), dSf_Haq(0:1), d2Saf_Haq(0:1), d2Sqf_Haq(0:1) ! Values at Half a & Half q (1D in f)
      real(wp) V(3,0:1,0:1,0:1),dV(3,3,0:1,0:1,0:1),d23V(4,3,0:1,0:1,0:1) ! Values at corners (3D in a,q,f)
      real(wp) V_Ha(3,0:1,0:1), dVq_Ha(3,0:1,0:1), dVf_Ha(3,0:1,0:1), d2Vqf_Ha(3,0:1,0:1)  ! Values at Half a (2D in q,f)
      real(wp) V_Haq(3,0:1), dVf_Haq(3,0:1) ! Values at Half a & Half q (1D in f)
      real(wp) F(0:1,0:1),d12F(3,0:1,0:1) ! Values at corners (2D in a,q)
      real(wp) F_Ha(0:1),dFa_Ha(0:1),dFq_Ha(0:1),dFaq_Ha(0:1) ! Values at Half a (1D in q)

      real(wp) in_a,in_gq,in_gf
      real(wp) P_R,P_grad_R(2),P_Z,P_grad_Z(2)
      real(wp) P_BB, P_grad_BB(3), P_B_cntv(3), P_B_cov(3), P_j_cntv(3)



      aid = max(1,min(dim_a-1, floor((in_a-a_1D(1))*inv_da +1)))
      qid = max(1,min(dim_gq-1, floor((in_gq-gq_1D(1))*inv_dgq +1)))
      fid = max(1,min(dim_gf-1, floor((in_gf-gf_1D(1))*inv_dgf +1)))

      at = (in_a-a_1D(aid))*inv_da
      qt = (in_gq-gq_1D(qid))*inv_dgq
      ft = (in_gf-gf_1D(fid))*inv_dgf


      tt = [at,qt,ft]
      du = [da,dgq,dgf]

      h00 = 2*tt**3 -3*tt**2 +1
      h10 = (tt**3 -2*tt**2 +tt)*du
      h01 = -2*tt**3 +3*tt**2
      h11 = (tt**3 -tt**2)*du
      dh00 = (6*tt**2 -6*tt)/du
      dh10 = 3*tt**2 -4*tt +1
      dh01 = (-6*tt**2 +6*tt)/du
      dh11 = 3*tt**2 -2*tt



      F  = R_2D(qid:qid+1,aid:aid+1)
      d12F = d12_R_2D(:,qid:qid+1,aid:aid+1)
      ! Find values at Half a (for 1D in gq)
      F_Ha = h00(1)*F(:,0) +h10(1)*d12F(1,:,0) +h01(1)*F(:,1) +h11(1)*d12F(1,:,1)
      dFa_Ha = dh00(1)*F(:,0) +dh10(1)*d12F(1,:,0) +dh01(1)*F(:,1) +dh11(1)*d12F(1,:,1)
      dFq_Ha = h00(1)*d12F(2,:,0) +h10(1)*d12F(3,:,0) +h01(1)*d12F(2,:,1) +h11(1)*d12F(3,:,1)
      dFaq_Ha = dh00(1)*d12F(2,:,0) +dh10(1)*d12F(3,:,0) +dh01(1)*d12F(2,:,1) +dh11(1)*d12F(3,:,1)
      ! Find values at Half a & Half gq (at the given point)
      P_R = h00(2)*F_Ha(0) +h10(2)*dFq_Ha(0) +h01(2)*F_Ha(1) +h11(2)*dFq_Ha(1) ! F
      P_grad_R(1) = h00(2)*dFa_Ha(0) +h10(2)*dFaq_Ha(0) +h01(2)*dFa_Ha(1) +h11(2)*dFaq_Ha(1) ! dFda
      P_grad_R(2) = dh00(2)*F_Ha(0) +dh10(2)*dFq_Ha(0) +dh01(2)*F_Ha(1) +dh11(2)*dFq_Ha(1) ! dFdgq


      F  = Z_2D(qid:qid+1,aid:aid+1)
      d12F = d12_Z_2D(:,qid:qid+1,aid:aid+1)
      ! Find values at Half a (for 1D in gq)
      F_Ha = h00(1)*F(:,0) +h10(1)*d12F(1,:,0) +h01(1)*F(:,1) +h11(1)*d12F(1,:,1)
      dFa_Ha = dh00(1)*F(:,0) +dh10(1)*d12F(1,:,0) +dh01(1)*F(:,1) +dh11(1)*d12F(1,:,1)
      dFq_Ha = h00(1)*d12F(2,:,0) +h10(1)*d12F(3,:,0) +h01(1)*d12F(2,:,1) +h11(1)*d12F(3,:,1)
      dFaq_Ha = dh00(1)*d12F(2,:,0) +dh10(1)*d12F(3,:,0) +dh01(1)*d12F(2,:,1) +dh11(1)*d12F(3,:,1)
      ! Find values at Half a & Half gq (at the given point)
      P_Z = h00(2)*F_Ha(0) +h10(2)*dFq_Ha(0) +h01(2)*F_Ha(1) +h11(2)*dFq_Ha(1) ! F
      P_grad_Z(1) = h00(2)*dFa_Ha(0) +h10(2)*dFaq_Ha(0) +h01(2)*dFa_Ha(1) +h11(2)*dFaq_Ha(1) ! dFda
      P_grad_Z(2) = dh00(2)*F_Ha(0) +dh10(2)*dFq_Ha(0) +dh01(2)*F_Ha(1) +dh11(2)*dFq_Ha(1) ! dFdgq




      ! Some Scalars
      S  = BB(qid:qid+1,aid:aid+1,fid:fid+1)
      dS = grad_BB(:,qid:qid+1,aid:aid+1,fid:fid+1)
      d23S = d23_BB(:,qid:qid+1,aid:aid+1,fid:fid+1)

      ! Find values at Half a (for 2D in gq,gf)
      S_Ha = h00(1)*S(:,0,:) +h10(1)*dS(1,:,0,:) +h01(1)*S(:,1,:) +h11(1)*dS(1,:,1,:)
      dSa_Ha=dh00(1)*S(:,0,:) +dh10(1)*dS(1,:,0,:) +dh01(1)*S(:,1,:) +dh11(1)*dS(1,:,1,:)
      dSq_Ha = h00(1)*dS(2,:,0,:) +h10(1)*d23S(1,:,0,:) +h01(1)*dS(2,:,1,:) +h11(1)*d23S(1,:,1,:)
      dSf_Ha = h00(1)*dS(3,:,0,:) +h10(1)*d23S(2,:,0,:) +h01(1)*dS(3,:,1,:) +h11(1)*d23S(2,:,1,:)
      d2Saq_Ha = dh00(1)*dS(2,:,0,:) +dh10(1)*d23S(1,:,0,:) +dh01(1)*dS(2,:,1,:) +dh11(1)*d23S(1,:,1,:)
      d2Saf_Ha = dh00(1)*dS(3,:,0,:) +dh10(1)*d23S(2,:,0,:) +dh01(1)*dS(3,:,1,:) +dh11(1)*d23S(2,:,1,:)
      d2Sqf_Ha = h00(1)*d23S(3,:,0,:) +h10(1)*d23S(4,:,0,:) +h01(1)*d23S(3,:,1,:) +h11(1)*d23S(4,:,1,:)
      d3Saqf_Ha=dh00(1)*d23S(3,:,0,:) +dh10(1)*d23S(4,:,0,:) +dh01(1)*d23S(3,:,1,:) +dh11(1)*d23S(4,:,1,:)
      ! Find values at Half a & Half gq (for 1D in gf)
      S_Haq = h00(2)*S_Ha(0,:) +h10(2)*dSq_Ha(0,:) +h01(2)*S_Ha(1,:) +h11(2)*dSq_Ha(1,:)
      dSa_Haq =  h00(2)*dSa_Ha(0,:) +h10(2)*d2Saq_Ha(0,:) +h01(2)*dSa_Ha(0,:) +h11(2)*d2Saq_Ha(0,:)
      dSq_Haq = dh00(2)*S_Ha(0,:) +dh10(2)*dSq_Ha(0,:) +dh01(2)*S_Ha(1,:) +dh11(2)*dSq_Ha(1,:)
      dSf_Haq =  h00(2)*dSf_Ha(0,:) +h10(2)*d2Sqf_Ha(0,:) +h01(2)*dSf_Ha(1,:) +h11(2)*d2Sqf_Ha(1,:)
      d2Saf_Haq = h00(2)*d2Saf_Ha(0,:) +h10(2)*d3Saqf_Ha(0,:) +h01(2)*d2Saf_Ha(1,:) +h11(2)*d3Saqf_Ha(1,:)
      d2Sqf_Haq = dh00(2)*dSf_Ha(0,:) +dh10(2)*d2Sqf_Ha(0,:) +dh01(2)*dSf_Ha(1,:) +dh11(2)*d2Sqf_Ha(1,:)
      ! Find interpolated value at the given point
      P_BB = h00(3)*S_Haq(0) +h10(3)*dSf_Haq(0) +h01(3)*S_Haq(1) +h11(3)*dSf_Haq(1)
      P_grad_BB(1) = h00(3)*dSa_Haq(0) +h10(3)*d2Saf_Haq(0) +h01(3)*dSa_Haq(1) +h11(3)*d2Saf_Haq(1)
      P_grad_BB(2) = h00(3)*dSq_Haq(0) +h10(3)*d2Sqf_Haq(0) +h01(3)*dSq_Haq(1) +h11(3)*d2Sqf_Haq(1)


      V  = B_cov(:,qid:qid+1,aid:aid+1,fid:fid+1)
      dV = grad_B_cov(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
      d23V = d23_B_cov(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
      ! Find values at Half a (for 2D in gq,gf)
      V_Ha = h00(1)*V(:,:,0,:) +h10(1)*dV(1,:,:,0,:) +h01(1)*V(:,:,1,:) +h11(1)*dV(1,:,:,1,:)
      dVq_Ha = h00(1)*dV(2,:,:,0,:) +h10(1)*d23V(1,:,:,0,:) +h01(1)*dV(2,:,:,1,:) +h11(1)*d23V(1,:,:,1,:)
      dVf_Ha = h00(1)*dV(3,:,:,0,:) +h10(1)*d23V(2,:,:,0,:) +h01(1)*dV(3,:,:,1,:) +h11(1)*d23V(2,:,:,1,:)
      d2Vqf_Ha = h00(1)*d23V(3,:,:,0,:) +h10(1)*d23V(4,:,:,0,:) +h01(1)*d23V(3,:,:,1,:) +h11(1)*d23V(4,:,:,1,:)
      ! Find values at Half a & Half gq (for 1D in gf)
      V_Haq = h00(2)*V_Ha(:,0,:) +h10(2)*dVq_Ha(:,0,:) +h01(2)*V_Ha(:,1,:) +h11(2)*dVq_Ha(:,1,:)
      dVf_Haq =  h00(2)*dVf_Ha(:,0,:) +h10(2)*d2Vqf_Ha(:,0,:) +h01(2)*dVf_Ha(:,1,:) +h11(2)*d2Vqf_Ha(:,1,:)
      ! Find interpolated value at the given point
      P_B_cov = h00(3)*V_Haq(:,0) +h10(3)*dVf_Haq(:,0) +h01(3)*V_Haq(:,1) +h11(3)*dVf_Haq(:,1)


      V  = B_cntv(:,qid:qid+1,aid:aid+1,fid:fid+1)
      dV = grad_B_cntv(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
      d23V = d23_B_cntv(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
      ! Find values at Half a (for 2D in gq,gf)
      V_Ha = h00(1)*V(:,:,0,:) +h10(1)*dV(1,:,:,0,:) +h01(1)*V(:,:,1,:) +h11(1)*dV(1,:,:,1,:)
      dVq_Ha = h00(1)*dV(2,:,:,0,:) +h10(1)*d23V(1,:,:,0,:) +h01(1)*dV(2,:,:,1,:) +h11(1)*d23V(1,:,:,1,:)
      dVf_Ha = h00(1)*dV(3,:,:,0,:) +h10(1)*d23V(2,:,:,0,:) +h01(1)*dV(3,:,:,1,:) +h11(1)*d23V(2,:,:,1,:)
      d2Vqf_Ha = h00(1)*d23V(3,:,:,0,:) +h10(1)*d23V(4,:,:,0,:) +h01(1)*d23V(3,:,:,1,:) +h11(1)*d23V(4,:,:,1,:)
      ! Find values at Half a & Half gq (for 1D in gf)
      V_Haq = h00(2)*V_Ha(:,0,:) +h10(2)*dVq_Ha(:,0,:) +h01(2)*V_Ha(:,1,:) +h11(2)*dVq_Ha(:,1,:)
      dVf_Haq =  h00(2)*dVf_Ha(:,0,:) +h10(2)*d2Vqf_Ha(:,0,:) +h01(2)*dVf_Ha(:,1,:) +h11(2)*d2Vqf_Ha(:,1,:)
      ! Find interpolated value at the given point
      P_B_cntv = h00(3)*V_Haq(:,0) +h10(3)*dVf_Haq(:,0) +h01(3)*V_Haq(:,1) +h11(3)*dVf_Haq(:,1)


      V  = j_cntv(:,qid:qid+1,aid:aid+1,fid:fid+1)
      dV = grad_j_cntv(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
      d23V = d23_j_cntv(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
      ! Find values at Half a (for 2D in gq,gf)
      V_Ha = h00(1)*V(:,:,0,:) +h10(1)*dV(1,:,:,0,:) +h01(1)*V(:,:,1,:) +h11(1)*dV(1,:,:,1,:)
      dVq_Ha = h00(1)*dV(2,:,:,0,:) +h10(1)*d23V(1,:,:,0,:) +h01(1)*dV(2,:,:,1,:) +h11(1)*d23V(1,:,:,1,:)
      dVf_Ha = h00(1)*dV(3,:,:,0,:) +h10(1)*d23V(2,:,:,0,:) +h01(1)*dV(3,:,:,1,:) +h11(1)*d23V(2,:,:,1,:)
      d2Vqf_Ha = h00(1)*d23V(3,:,:,0,:) +h10(1)*d23V(4,:,:,0,:) +h01(1)*d23V(3,:,:,1,:) +h11(1)*d23V(4,:,:,1,:)
      ! Find values at Half a & Half gq (for 1D in gf)
      V_Haq = h00(2)*V_Ha(:,0,:) +h10(2)*dVq_Ha(:,0,:) +h01(2)*V_Ha(:,1,:) +h11(2)*dVq_Ha(:,1,:)
      dVf_Haq =  h00(2)*dVf_Ha(:,0,:) +h10(2)*d2Vqf_Ha(:,0,:) +h01(2)*dVf_Ha(:,1,:) +h11(2)*d2Vqf_Ha(:,1,:)
      ! Find interpolated value at the given point
      P_j_cntv = h00(3)*V_Haq(:,0) +h10(3)*dVf_Haq(:,0) +h01(3)*V_Haq(:,1) +h11(3)*dVf_Haq(:,1)

   end subroutine Cubic_Interpolation_of_VariousStuff_for_GCmotion




   subroutine Cubic_Interpolation_of_VariousStuff_for_GCmotion_at_zeon(P_R,P_grad_R,P_Z,P_grad_Z,P_BB,P_grad_BB,P_B_cov,P_B_cntv,P_j_cntv)
      use global_parameters,only:me
      use particle_array
      implicit none
      ! real(wp),intent(in) :: in_a,in_gq,in_gf
      ! real(wp),intent(in) :: in_Vec(3,dim_gq,dim_a,dim_gf)
      ! real(wp),intent(in) :: in_gradVec(3,3,dim_gq,dim_a,dim_gf), in_d23Vec(4,3,dim_gq,dim_a,dim_gf)
      ! real(wp),intent(out) :: out_Vec(3)
      real(wp) h00(3),h10(3),h01(3),h11(3)
      real(wp) dh00(3),dh10(3),dh01(3),dh11(3)
      integer aid,qid,fid
      real(wp) at,qt,ft,tt(3),du(3)
      real(wp) S(0:1,0:1,0:1),dS(3,0:1,0:1,0:1),d23S(4,0:1,0:1,0:1) ! Values at corners (3D in a,q,f)
      real(wp) S_Ha(0:1,0:1), dSa_Ha(0:1,0:1), dSq_Ha(0:1,0:1), dSf_Ha(0:1,0:1), d2Saq_Ha(0:1,0:1), d2Saf_Ha(0:1,0:1), d2Sqf_Ha(0:1,0:1), d3Saqf_Ha(0:1,0:1)  ! Values at Half a (2D in q,f)
      real(wp) S_Haq(0:1), dSa_Haq(0:1), dSq_Haq(0:1), dSf_Haq(0:1), d2Saf_Haq(0:1), d2Sqf_Haq(0:1) ! Values at Half a & Half q (1D in f)
      real(wp) V(3,0:1,0:1,0:1),dV(3,3,0:1,0:1,0:1),d23V(4,3,0:1,0:1,0:1) ! Values at corners (3D in a,q,f)
      real(wp) V_Ha(3,0:1,0:1), dVq_Ha(3,0:1,0:1), dVf_Ha(3,0:1,0:1), d2Vqf_Ha(3,0:1,0:1)  ! Values at Half a (2D in q,f)
      real(wp) V_Haq(3,0:1), dVf_Haq(3,0:1) ! Values at Half a & Half q (1D in f)
      real(wp) F(0:1,0:1),d12F(3,0:1,0:1) ! Values at corners (2D in a,q)
      real(wp) F_Ha(0:1),dFa_Ha(0:1),dFq_Ha(0:1),dFaq_Ha(0:1) ! Values at Half a (1D in q)

      real(wp) in_a,in_gq,in_gf
      real(wp) P_R(me),P_grad_R(2,me),P_Z(me),P_grad_Z(2,me)
      real(wp) P_BB(me), P_grad_BB(3,me), P_B_cntv(3,me), P_B_cov(3,me), P_j_cntv(3,me)

      integer m

      !$acc parallel loop gang vector present(zeon)
      do m=1,me

         in_a = zeon(1,m)
         in_gq = zeon(2,m)
         in_gf = zeon(3,m)

         aid = max(1,min(dim_a-1, floor((in_a-a_1D(1))*inv_da +1)))
         qid = max(1,min(dim_gq-1, floor((in_gq-gq_1D(1))*inv_dgq +1)))
         fid = max(1,min(dim_gf-1, floor((in_gf-gf_1D(1))*inv_dgf +1)))

         at = (in_a-a_1D(aid))*inv_da
         qt = (in_gq-gq_1D(qid))*inv_dgq
         ft = (in_gf-gf_1D(fid))*inv_dgf


         tt = [at,qt,ft]
         du = [da,dgq,dgf]

         h00 = 2*tt**3 -3*tt**2 +1
         h10 = (tt**3 -2*tt**2 +tt)*du
         h01 = -2*tt**3 +3*tt**2
         h11 = (tt**3 -tt**2)*du
         dh00 = (6*tt**2 -6*tt)/du
         dh10 = 3*tt**2 -4*tt +1
         dh01 = (-6*tt**2 +6*tt)/du
         dh11 = 3*tt**2 -2*tt


         F  = R_2D(qid:qid+1,aid:aid+1)
         d12F = d12_R_2D(:,qid:qid+1,aid:aid+1)
         ! Find values at Half a (for 1D in gq)
         F_Ha = h00(1)*F(:,0) +h10(1)*d12F(1,:,0) +h01(1)*F(:,1) +h11(1)*d12F(1,:,1)
         dFa_Ha = dh00(1)*F(:,0) +dh10(1)*d12F(1,:,0) +dh01(1)*F(:,1) +dh11(1)*d12F(1,:,1)
         dFq_Ha = h00(1)*d12F(2,:,0) +h10(1)*d12F(3,:,0) +h01(1)*d12F(2,:,1) +h11(1)*d12F(3,:,1)
         dFaq_Ha = dh00(1)*d12F(2,:,0) +dh10(1)*d12F(3,:,0) +dh01(1)*d12F(2,:,1) +dh11(1)*d12F(3,:,1)
         ! Find values at Half a & Half gq (at the given point)
         P_R(m) = h00(2)*F_Ha(0) +h10(2)*dFq_Ha(0) +h01(2)*F_Ha(1) +h11(2)*dFq_Ha(1) ! F
         P_grad_R(1,m) = h00(2)*dFa_Ha(0) +h10(2)*dFaq_Ha(0) +h01(2)*dFa_Ha(1) +h11(2)*dFaq_Ha(1) ! dFda
         P_grad_R(2,m) = dh00(2)*F_Ha(0) +dh10(2)*dFq_Ha(0) +dh01(2)*F_Ha(1) +dh11(2)*dFq_Ha(1) ! dFdgq


         F  = Z_2D(qid:qid+1,aid:aid+1)
         d12F = d12_Z_2D(:,qid:qid+1,aid:aid+1)
         ! Find values at Half a (for 1D in gq)
         F_Ha = h00(1)*F(:,0) +h10(1)*d12F(1,:,0) +h01(1)*F(:,1) +h11(1)*d12F(1,:,1)
         dFa_Ha = dh00(1)*F(:,0) +dh10(1)*d12F(1,:,0) +dh01(1)*F(:,1) +dh11(1)*d12F(1,:,1)
         dFq_Ha = h00(1)*d12F(2,:,0) +h10(1)*d12F(3,:,0) +h01(1)*d12F(2,:,1) +h11(1)*d12F(3,:,1)
         dFaq_Ha = dh00(1)*d12F(2,:,0) +dh10(1)*d12F(3,:,0) +dh01(1)*d12F(2,:,1) +dh11(1)*d12F(3,:,1)
         ! Find values at Half a & Half gq (at the given point)
         P_Z(m) = h00(2)*F_Ha(0) +h10(2)*dFq_Ha(0) +h01(2)*F_Ha(1) +h11(2)*dFq_Ha(1) ! F
         P_grad_Z(1,m) = h00(2)*dFa_Ha(0) +h10(2)*dFaq_Ha(0) +h01(2)*dFa_Ha(1) +h11(2)*dFaq_Ha(1) ! dFda
         P_grad_Z(2,m) = dh00(2)*F_Ha(0) +dh10(2)*dFq_Ha(0) +dh01(2)*F_Ha(1) +dh11(2)*dFq_Ha(1) ! dFdgq




         ! Some Scalars
         S  = BB(qid:qid+1,aid:aid+1,fid:fid+1)
         dS = grad_BB(:,qid:qid+1,aid:aid+1,fid:fid+1)
         d23S = d23_BB(:,qid:qid+1,aid:aid+1,fid:fid+1)

         ! Find values at Half a (for 2D in gq,gf)
         S_Ha = h00(1)*S(:,0,:) +h10(1)*dS(1,:,0,:) +h01(1)*S(:,1,:) +h11(1)*dS(1,:,1,:)
         dSa_Ha=dh00(1)*S(:,0,:) +dh10(1)*dS(1,:,0,:) +dh01(1)*S(:,1,:) +dh11(1)*dS(1,:,1,:)
         dSq_Ha = h00(1)*dS(2,:,0,:) +h10(1)*d23S(1,:,0,:) +h01(1)*dS(2,:,1,:) +h11(1)*d23S(1,:,1,:)
         dSf_Ha = h00(1)*dS(3,:,0,:) +h10(1)*d23S(2,:,0,:) +h01(1)*dS(3,:,1,:) +h11(1)*d23S(2,:,1,:)
         d2Saq_Ha = dh00(1)*dS(2,:,0,:) +dh10(1)*d23S(1,:,0,:) +dh01(1)*dS(2,:,1,:) +dh11(1)*d23S(1,:,1,:)
         d2Saf_Ha = dh00(1)*dS(3,:,0,:) +dh10(1)*d23S(2,:,0,:) +dh01(1)*dS(3,:,1,:) +dh11(1)*d23S(2,:,1,:)
         d2Sqf_Ha = h00(1)*d23S(3,:,0,:) +h10(1)*d23S(4,:,0,:) +h01(1)*d23S(3,:,1,:) +h11(1)*d23S(4,:,1,:)
         d3Saqf_Ha=dh00(1)*d23S(3,:,0,:) +dh10(1)*d23S(4,:,0,:) +dh01(1)*d23S(3,:,1,:) +dh11(1)*d23S(4,:,1,:)
         ! Find values at Half a & Half gq (for 1D in gf)
         S_Haq = h00(2)*S_Ha(0,:) +h10(2)*dSq_Ha(0,:) +h01(2)*S_Ha(1,:) +h11(2)*dSq_Ha(1,:)
         dSa_Haq =  h00(2)*dSa_Ha(0,:) +h10(2)*d2Saq_Ha(0,:) +h01(2)*dSa_Ha(0,:) +h11(2)*d2Saq_Ha(0,:)
         dSq_Haq = dh00(2)*S_Ha(0,:) +dh10(2)*dSq_Ha(0,:) +dh01(2)*S_Ha(1,:) +dh11(2)*dSq_Ha(1,:)
         dSf_Haq =  h00(2)*dSf_Ha(0,:) +h10(2)*d2Sqf_Ha(0,:) +h01(2)*dSf_Ha(1,:) +h11(2)*d2Sqf_Ha(1,:)
         d2Saf_Haq = h00(2)*d2Saf_Ha(0,:) +h10(2)*d3Saqf_Ha(0,:) +h01(2)*d2Saf_Ha(1,:) +h11(2)*d3Saqf_Ha(1,:)
         d2Sqf_Haq = dh00(2)*dSf_Ha(0,:) +dh10(2)*d2Sqf_Ha(0,:) +dh01(2)*dSf_Ha(1,:) +dh11(2)*d2Sqf_Ha(1,:)
         ! Find interpolated value at the given point
         P_BB(m) = h00(3)*S_Haq(0) +h10(3)*dSf_Haq(0) +h01(3)*S_Haq(1) +h11(3)*dSf_Haq(1)
         P_grad_BB(1,m) = h00(3)*dSa_Haq(0) +h10(3)*d2Saf_Haq(0) +h01(3)*dSa_Haq(1) +h11(3)*d2Saf_Haq(1)
         P_grad_BB(2,m) = h00(3)*dSq_Haq(0) +h10(3)*d2Sqf_Haq(0) +h01(3)*dSq_Haq(1) +h11(3)*d2Sqf_Haq(1)


         V  = B_cov(:,qid:qid+1,aid:aid+1,fid:fid+1)
         dV = grad_B_cov(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
         d23V = d23_B_cov(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
         ! Find values at Half a (for 2D in gq,gf)
         V_Ha = h00(1)*V(:,:,0,:) +h10(1)*dV(1,:,:,0,:) +h01(1)*V(:,:,1,:) +h11(1)*dV(1,:,:,1,:)
         dVq_Ha = h00(1)*dV(2,:,:,0,:) +h10(1)*d23V(1,:,:,0,:) +h01(1)*dV(2,:,:,1,:) +h11(1)*d23V(1,:,:,1,:)
         dVf_Ha = h00(1)*dV(3,:,:,0,:) +h10(1)*d23V(2,:,:,0,:) +h01(1)*dV(3,:,:,1,:) +h11(1)*d23V(2,:,:,1,:)
         d2Vqf_Ha = h00(1)*d23V(3,:,:,0,:) +h10(1)*d23V(4,:,:,0,:) +h01(1)*d23V(3,:,:,1,:) +h11(1)*d23V(4,:,:,1,:)
         ! Find values at Half a & Half gq (for 1D in gf)
         V_Haq = h00(2)*V_Ha(:,0,:) +h10(2)*dVq_Ha(:,0,:) +h01(2)*V_Ha(:,1,:) +h11(2)*dVq_Ha(:,1,:)
         dVf_Haq =  h00(2)*dVf_Ha(:,0,:) +h10(2)*d2Vqf_Ha(:,0,:) +h01(2)*dVf_Ha(:,1,:) +h11(2)*d2Vqf_Ha(:,1,:)
         ! Find interpolated value at the given point
         P_B_cov(:,m) = h00(3)*V_Haq(:,0) +h10(3)*dVf_Haq(:,0) +h01(3)*V_Haq(:,1) +h11(3)*dVf_Haq(:,1)


         V  = B_cntv(:,qid:qid+1,aid:aid+1,fid:fid+1)
         dV = grad_B_cntv(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
         d23V = d23_B_cntv(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
         ! Find values at Half a (for 2D in gq,gf)
         V_Ha = h00(1)*V(:,:,0,:) +h10(1)*dV(1,:,:,0,:) +h01(1)*V(:,:,1,:) +h11(1)*dV(1,:,:,1,:)
         dVq_Ha = h00(1)*dV(2,:,:,0,:) +h10(1)*d23V(1,:,:,0,:) +h01(1)*dV(2,:,:,1,:) +h11(1)*d23V(1,:,:,1,:)
         dVf_Ha = h00(1)*dV(3,:,:,0,:) +h10(1)*d23V(2,:,:,0,:) +h01(1)*dV(3,:,:,1,:) +h11(1)*d23V(2,:,:,1,:)
         d2Vqf_Ha = h00(1)*d23V(3,:,:,0,:) +h10(1)*d23V(4,:,:,0,:) +h01(1)*d23V(3,:,:,1,:) +h11(1)*d23V(4,:,:,1,:)
         ! Find values at Half a & Half gq (for 1D in gf)
         V_Haq = h00(2)*V_Ha(:,0,:) +h10(2)*dVq_Ha(:,0,:) +h01(2)*V_Ha(:,1,:) +h11(2)*dVq_Ha(:,1,:)
         dVf_Haq =  h00(2)*dVf_Ha(:,0,:) +h10(2)*d2Vqf_Ha(:,0,:) +h01(2)*dVf_Ha(:,1,:) +h11(2)*d2Vqf_Ha(:,1,:)
         ! Find interpolated value at the given point
         P_B_cntv(:,m) = h00(3)*V_Haq(:,0) +h10(3)*dVf_Haq(:,0) +h01(3)*V_Haq(:,1) +h11(3)*dVf_Haq(:,1)


         V  = j_cntv(:,qid:qid+1,aid:aid+1,fid:fid+1)
         dV = grad_j_cntv(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
         d23V = d23_j_cntv(:,:,qid:qid+1,aid:aid+1,fid:fid+1)
         ! Find values at Half a (for 2D in gq,gf)
         V_Ha = h00(1)*V(:,:,0,:) +h10(1)*dV(1,:,:,0,:) +h01(1)*V(:,:,1,:) +h11(1)*dV(1,:,:,1,:)
         dVq_Ha = h00(1)*dV(2,:,:,0,:) +h10(1)*d23V(1,:,:,0,:) +h01(1)*dV(2,:,:,1,:) +h11(1)*d23V(1,:,:,1,:)
         dVf_Ha = h00(1)*dV(3,:,:,0,:) +h10(1)*d23V(2,:,:,0,:) +h01(1)*dV(3,:,:,1,:) +h11(1)*d23V(2,:,:,1,:)
         d2Vqf_Ha = h00(1)*d23V(3,:,:,0,:) +h10(1)*d23V(4,:,:,0,:) +h01(1)*d23V(3,:,:,1,:) +h11(1)*d23V(4,:,:,1,:)
         ! Find values at Half a & Half gq (for 1D in gf)
         V_Haq = h00(2)*V_Ha(:,0,:) +h10(2)*dVq_Ha(:,0,:) +h01(2)*V_Ha(:,1,:) +h11(2)*dVq_Ha(:,1,:)
         dVf_Haq =  h00(2)*dVf_Ha(:,0,:) +h10(2)*d2Vqf_Ha(:,0,:) +h01(2)*dVf_Ha(:,1,:) +h11(2)*d2Vqf_Ha(:,1,:)
         ! Find interpolated value at the given point
         P_j_cntv(:,m) = h00(3)*V_Haq(:,0) +h10(3)*dVf_Haq(:,0) +h01(3)*V_Haq(:,1) +h11(3)*dVf_Haq(:,1)

      enddo
   end subroutine Cubic_Interpolation_of_VariousStuff_for_GCmotion_at_zeon









end module Magnetic_Grid