module MGY_PETSc_Interface_Data
#include <petsc/finclude/petsc.h>
   use petscksp
   use precision, only:wp
   implicit none
   type PEST_Indices
      PetscInt gqid, aid, gfid
   end type PEST_Indices

   ! Boundary Conditions
   ! (BC0 and BC1 are inner (s=s0) and outer (s=s1) conditions, respectively)
   PetscInt :: BC0_type=0, BC1_type=0
   PetscReal :: BC0_value=0,  BC1_value=0

   type(PEST_Indices),allocatable,dimension(:) :: mapG2P ! From Global node index to PEST node indices
   PetscInt, allocatable, dimension(:,:,:) :: mapP2G, mapP2L ! From PEST indices to Global or Local node index
   PetscInt, allocatable, dimension(:) :: mapL2G ! From Local indices to Global indices (this depends on the MPI process)
   PetscInt Lsize, Gsize, sum_Lsize_partd ! Size of Local nodes and Global nodes
   PetscInt, allocatable, dimension(:) :: recv_Lsize_partd, displs_Lsize_partd
   PetscInt, allocatable, dimension(:) :: recv_mapL2G_partd
   PetscInt Gid_start, Gid_end
   Mat PMat ! GK Poisson equation's LHS matrix
   Mat AMat ! GK Ampere equation's LHS matrix
   Vec bVec, xVec, ResVec
   KSP ksp_P ! ksp for GK Poisson equation
   KSP ksp_A ! ksp for GK Ampere equation

   ! Follolwing (Jmat, yVec) are for solving nonlinear Poisson equation in the case of Flag_dne_adiab==3
   !Jacobian Matrix (PMat + (ne0*e^2/Te0)*exp(e*phi/Te0)), in the case of Flag_dne_adiab==3 (nonlinear adiabatic response)
   Mat Jmat, pre_Jmat !pre_Jamt is a preconditioner of Jmat which will be updated occasionaly to reduce the computational time
   ! F vector (a system of equations) for Newton-Raphson method
   Vec yVec ! yVec = (Jmat)^(-1)*(FVec)
   Vec srcVec ! source (RHS) of the nonlinear Poisson equation (e*(ni0+dni-dneh)) to calculate residual

   PetscInt, allocatable, dimension(:,:) :: tarGid
   PetscReal, allocatable, dimension(:,:) :: tarVal

   integer, allocatable :: o_nnz(:), d_nnz(:) !o_nnz and d_nnz are number of off-diagonal or diagonal nonzero terms
   integer max_nnz ! maximum number of nonzero compoentns per each row

   ! Following Tensors will be used to determine PMat components (Poisson equation)
   ! Note that they are composed with the contravariant components
   ! Metric_T:  contravariant Metrics; g^(ij)=e^i*e^j, where e is covariant vector on the PEST coordinate
   ! Polar_T: Polarization shielding tensor. (n_s*m_s/B^2)*(I_tensor-b*b), where b is unit vector of the magnetic field
   ! Eps_T: Metric_T + Polar_T
   ! CT: Jacobian*Eps_T/(dui*duj)

   ! Metric_T is static once it is set on a given PEST structred coordinate
   ! Metric_T is used for solving both GK Poisson and GK Ampere equation
   ! Metric_T(i,j) = g^ij 
   real(wp), allocatable, dimension(:,:,:,:) :: Metric_T! On the poloidal plane (PEST structured coordinate)

   ! Polar_T and Eps_T are used only for the GK Poisson equation
   ! Polar_T = Sigma_over_s(ns*ms/B^2 * (Metric_T - B^i*B^j/B^2))
   ! Eps_T=eps0*Metric_T+Polar_T (Flag=0) or Eps_T=eps0*Metric_T (Flag=1) or Eps_T=Polar_T (Flag=2)
   real(wp), allocatable :: Polar_T(:,:,:,:), Eps_T(:,:,:,:,:) 


   ! CT is Coefficient Tensor that is needed to solve laplacian-like operator (nabla_perp*nabla_perp) in GK Poisson or GK Ampere equations 
   ! For GK Poisson equation, CT = J*Eps_T(i,j)/(du^i*du^j), where J is jacobian
   ! For GK Ampere equation, CT= J*(Metric_T(i,j)-B^i*B^j/B^2)/(du^i*du^j)
   ! CT values are updated at every time when solving either GK Poisson or GK Ampere equation
   real(wp), allocatable, dimension(:,:,:,:,:) :: CT ! LHS coefficient tensor including neighboring left and right poloidal planes


contains

subroutine Initialize_MGY_PETSc_Interface_Data(dim_gq,dim_a,dim_gf)
   use precision
   use global_parameters,only:stdout,mype,numberpe
   use particle_decomp
   implicit none
   integer, intent(in):: dim_gq,dim_a,dim_gf
   integer :: gqid, aid, gfid, Gid, Lid, tnum, extra, i, j
   PetscInt neq_Pol2D
   real(wp) :: tt0,tt1
   PetscErrorCode ierr
   PetscScalar, pointer :: ptr_bVec(:)

   if (mype==0) write(stdout,*) "\n\n\n\n\n"
   if (mype==0) write(stdout,*) "=============================================="
   if (mype==0) write(stdout,*) "   START! Initialize_MGY_PETSc_Interface_Data"
   if (mype==0) write(stdout,*) "=============================================="
   if (mype==0) write(stdout,*) "\n\n"
   call MPI_Barrier(PETSC_COMM_WORLD,ierr)



   Gsize = (dim_gq-1)*dim_a*dim_gf

   ! Manual partitioning for each poloidal plane
   neq_Pol2D = (dim_gq-1)*dim_a
   Lsize = neq_Pol2D/npartdom
   extra = neq_Pol2D-Lsize*npartdom

   Gid_start = myrank_toroidal*neq_Pol2D + Lsize*myrank_partd + min(myrank_partd,extra)
   Gid_end   = myrank_toroidal*neq_Pol2D + Lsize*(myrank_partd+1) + min(myrank_partd+1,extra)-1
   Lsize = Gid_end-Gid_start+1


   allocate(mapL2G(Lsize))
   allocate(mapG2P(0:Gsize-1)) ! 1-dimensional array (starting index is 0 because of PETSc's C-style indexing)
   allocate(mapP2G(0:dim_gq,0:dim_a+1,0:dim_gf+1)) ! 3-dimensional array (Note that 0, and dim+1 are ghost nodes)
   allocate(mapP2L(0:dim_gq,0:dim_a+1,0:dim_gf+1))

   ! Construct mapping from Local to Global index
   do Lid=1,Lsize
      mapL2G(Lid)=Gid_start+(Lid-1)
   enddo

   Gid=0 ! Starting index for global node is 0 because of PETSc's C-style indexing
   do gfid=1,dim_gf
      do aid=1,dim_a
         do gqid=1,dim_gq-1 !Last gqid is not used since it is the same as gqid=0 because of periodicity
            mapG2P(Gid)%gqid=gqid
            mapG2P(Gid)%aid=aid
            mapG2P(Gid)%gfid=gfid
            mapP2G(gqid,aid,gfid)=Gid
            Gid=Gid+1
         enddo
      enddo
   enddo

   mapP2L=-100 ! Default value is negative (-100), so one can tell whether a specific PETS point is out of the num_its process (mapP2L<0)
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      mapP2L(mapG2P(Gid)%gqid,mapG2P(Gid)%aid,mapG2P(Gid)%gfid)=Lid
   enddo

   ! Fill ghost nodes of mapP2G
   mapP2G(dim_gq,:,:)   =mapP2G(1,:,:)
   mapP2G(0,:,:)        =mapP2G(dim_gq-1,:,:)
   mapP2G(:,0,:)        =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   mapP2G(:,dim_a+1,:)  =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   mapP2G(:,:,0)        =mapP2G(:,:,dim_gf)
   mapP2G(:,:,dim_gf+1) =mapP2G(:,:,1)


   ! Construct recv_Lsize_partd which is the Lsize of the processors in each poloidal plane
   allocate(recv_Lsize_partd(npartdom), displs_Lsize_partd(npartdom))
   call MPI_ALLGATHER(Lsize,1,MPI_INTEGER,recv_Lsize_partd,1,MPI_INTEGER,partd_comm,ierr)
   displs_Lsize_partd(1)=0
   do i=1,npartdom-1
      displs_Lsize_partd(i+1)=displs_Lsize_partd(i)+recv_Lsize_partd(i)
   enddo
   call MPI_ALLGATHER(Lsize,1,MPI_INTEGER,recv_Lsize_partd,1,MPI_INTEGER,partd_comm,ierr)
   sum_Lsize_partd = sum(recv_Lsize_partd)
   allocate(recv_mapL2G_partd(sum_Lsize_partd))
   call MPI_ALLGATHERV(mapL2G,Lsize,MPI_INTEGER,recv_mapL2G_partd,recv_Lsize_partd,displs_Lsize_partd,MPI_INTEGER,partd_comm,ierr)


   if(mype==0) then
      write(stdout,*) "recv_Lsize_partd=",recv_Lsize_partd
      write(stdout,*) "displs_Lsize_partd=",displs_Lsize_partd
      call flush(stdout)
   endif

   ! call VecSet(bVec,5.0,ierr); CHKERRQ(ierr)
   ! call VecGetArrayF90(bVec,ptr_bVec,ierr)
   ! ptr_bVec(1) = 123.0
   ! call VecRestoreArrayF90(bVec,ptr_bVec,ierr)
   ! call VecView(bVec,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   write(stdout,"(3(A,I8))") "mype=",mype,", Gid_start=",Gid_start,", Gid_end=",Gid_end
   call flush(stdout)

   ! write(stdout,"(2(A,I8),1x,I8,A,I8)") "Gsize=",Gsize,", Lsize=",Lsize, 7*Lsize, ", sum(d_nnz+o_nnz)=",sum(d_nnz+o_nnz)
   ! call flush(stdout)
   ! call MatView(PMat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

   if(mype==0) then
      open(unit=200205,file="mapG2P.dat",action='write')
      do Gid=0,Gsize-1
         write(200205,"(4(I8))") Gid, mapG2P(Gid)%gqid, mapG2P(Gid)%aid, mapG2P(Gid)%gfid
      enddo
      call flush(200205)
      close(200205)

      open(unit=200205,file="mapL2G.dat",action='write')
      do Lid=1,Lsize
         write(200205,"(2(I8))") Lid, mapL2G(Lid)
      enddo
      call flush(200205)
      close(200205)


      open(unit=200205,file="mapP2LandP2G.dat",action='write')
      do gfid=0,dim_gf+1
         do aid=0,dim_a+1
            do gqid=0,dim_gq
               Gid=mapP2G(gqid,aid,gfid)
               write(200205,"(5(I8))") mapP2L(gqid,aid,gfid), Gid, mapG2P(Gid)%gqid,mapG2P(Gid)%aid,mapG2P(Gid)%gfid
            enddo
         enddo
      enddo
      call flush(200205)
      close(200205)

   endif


   call Construct_Metric_Tensor_on_Structured_Grid

end subroutine Initialize_MGY_PETSc_Interface_Data

subroutine Finalize_MGY_PETSc_Interface_Data
   implicit none

   deallocate(mapG2P,mapP2G,mapP2L,mapL2G)
   deallocate(recv_Lsize_partd,displs_Lsize_partd)
   deallocate(recv_mapL2G_partd)
   deallocate(d_nnz, o_nnz)

   if(allocated(Metric_T)) deallocate(Metric_T)
   if(allocated(Polar_T)) deallocate(Polar_T, Eps_T) ! for GK Poisson Equation
   if(allocated(CT)) deallocate(CT)

end subroutine Finalize_MGY_PETSc_Interface_Data


subroutine Construct_Metric_Tensor_on_Structured_Grid
   use precision, only:wp
   use Structured_Grid
   use particle_decomp
   implicit none

   integer sid,gqid,i,j,ierr

   if(.not.allocated(Metric_T)) allocate(Metric_T(3,3,0:SG%dim_gq,0:SG%dim_a+1))

   ! Metric Tensor is based on the static equilibrium magnetic field (PEST coordinate)
   do sid=1,SG%dim_a
      do gqid=1,SG%dim_gq
         Metric_T(:,:,gqid,sid)=0.0
         Metric_T(1,1,gqid,sid)=SG%g_ss_cntv(gqid,sid)
         Metric_T(2,1,gqid,sid)=SG%g_sq_cntv(gqid,sid)
         Metric_T(1,2,gqid,sid)=SG%g_sq_cntv(gqid,sid)
         Metric_T(2,2,gqid,sid)=SG%g_qq_cntv(gqid,sid)
         Metric_T(3,3,gqid,sid)=SG%g_ff_cntv(gqid,sid)
      enddo
   enddo

   ! Fill ghost nodes of Metric_T
   Metric_T(:,:,SG%dim_gq,:)=Metric_T(:,:,1,:)
   Metric_T(:,:,0,:)        =Metric_T(:,:,SG%dim_gq-1,:)
   Metric_T(:,:,:,0)        =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   Metric_T(:,:,:,SG%dim_a+1)  =-100 ! Assign Negative initial values which will be ignored in MatSetValues

end subroutine Construct_Metric_Tensor_on_Structured_Grid


subroutine Construct_SimplePoisson3D_Solver_using_Structured_Grid
   use Structured_Grid
   use particle_decomp
   use global_parameters,only:mype,stdout,zetamax,zetamin,mzetamax
   implicit none

   PetscInt Lid,Gid,gqid,aid,gfid,tnum
   PetscReal Coeff_gq, Coeff_a, Coeff_gf
   real tt0, tt1
   PetscErrorCode ierr
   real, external:: dsda

   ! call Initialize_MGY_PETSc_Interface_Data(SG%dim_gq,SG%dim_a,ntoroidal)

   max_nnz=7
   allocate(tarGid(max_nnz,Lsize), tarVal(max_nnz,Lsize))

   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   do Lid=1,Lsize
      Gid=mapL2G(Lid)

      gqid=mapG2P(Gid)%gqid
      aid =mapG2P(Gid)%aid
      gfid=mapG2P(Gid)%gfid

      Coeff_gq = SG%g_qq_cntv(gqid,aid)*SG%inv_dgq**2.0
      Coeff_a  = SG%g_aa_cntv(gqid,aid)*(dsda(SG%a_1D(aid))*SG%inv_ds)**2.0
      Coeff_gf = SG%g_ff_cntv(gqid,aid)/((zetamax-zetamin)**2.0)

      tnum=1
      tarGid(tnum,Lid)=mapP2G(gqid,aid,gfid);   tarVal(tnum,Lid)=2.0*(Coeff_gq+Coeff_a+Coeff_gf);  tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid-1,aid,gfid); tarVal(tnum,Lid)=-Coeff_gq; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid+1,aid,gfid); tarVal(tnum,Lid)=-Coeff_gq; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid,aid-1,gfid); tarVal(tnum,Lid)=-Coeff_a; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid,aid+1,gfid); tarVal(tnum,Lid)=-Coeff_a; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid,aid,gfid-1); tarVal(tnum,Lid)=-Coeff_gf; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid,aid,gfid+1); tarVal(tnum,Lid)=-Coeff_gf; tnum=tnum+1;
   enddo

   ! Apply Neumann (+Dirichlet) Condition at inner and outer boundary
   do gfid=1,ntoroidal
      do gqid=1,SG%dim_gq-1
         ! Neumann condition (Er=0 at inner, outer boundary)
         Lid=mapP2L(gqid,2,gfid)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=-1.0 !phi(aid=2)
            tarVal(4,Lid)=+1.0 !phi(aid=1)
         endif
         Lid=mapP2L(gqid,SG%dim_a-1,gfid)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=-1.0 !phi(aid=dim_a-1)
            tarVal(5,Lid)=+1.0 !phi(aid=dim_a)
         endif

         ! Neumann condition (E_theta=0 at inner, outer boundary)
         Lid=mapP2L(gqid,SG%dim_a,gfid)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=-1.0 !phi(aid=dim_a, gqid)
            tarVal(3,Lid)=+1.0 !phi(aid=dim_a, gqid+1)
         endif

         !Dirichlet condition (phi=0 at inner boundary)
         Lid=mapP2L(gqid,1,gfid)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=1.0 !phi(aid=1)
         endif
      enddo
   enddo

   ! To calculate numbers of non-zero elements for diagonal(d_nnz) & off-diagonal(o_nnz) blocks
   allocate(d_nnz(Lsize), o_nnz(Lsize))
   d_nnz=0; o_nnz=0;
   do Lid=1,Lsize
      do tnum=1,max_nnz
         if(tarGid(tnum,Lid)>=0) then
            if(tarGid(tnum,Lid)<Gid_start .or. tarGid(tnum,Lid)>Gid_end) then
               o_nnz(Lid)=o_nnz(Lid)+1
            else
               d_nnz(Lid)=d_nnz(Lid)+1
            endif
         endif
      enddo
   enddo

   ! Construct Matrix (PMat)
   call MatCreate(PETSC_COMM_WORLD,PMat,ierr); CHKERRQ(ierr)
   call MatSetSizes(PMat,Lsize,Lsize,Gsize,Gsize,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(PMat,ierr); CHKERRQ(ierr)
   ! Preallocate the memory using o_nnz and d_nnz
   call MatMPIAIJSetPreallocation(PMat,max_nnz,d_nnz,max_nnz,o_nnz,ierr); CHKERRQ(ierr)
   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(PMat,1,Gid,max_nnz,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo


   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call MatView(PMat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

   call VecCreate(PETSC_COMM_WORLD,bVec,ierr); CHKERRQ(ierr)
   call VecSetSizes(bVec,Lsize,Gsize,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(bVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,xVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,ResVec,ierr); CHKERRQ(ierr)


   ! Construct KSP
   tt0=MPI_WTIME()
   call KSPCreate(PETSC_COMM_WORLD,ksp_P,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp_P,PMat,PMat,ierr); CHKERRQ(ierr)
   ! call KSPSetTolerances(ksp_P,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp_P,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "ksp_P (GK Poisson) Create takes ",tt1-tt0," secs"
   call flush(stdout)


end subroutine Construct_SimplePoisson3D_Solver_using_Structured_Grid



subroutine Solve_Poisson3D_on_SG()
   use precision, only:wp
   ! use global_parameters, only:stdout,mype,ifelectron,mzetamax,mgrid
   use global_parameters
   use field_array

   use Structured_Grid
   use particle_decomp
   implicit none

   PetscScalar, pointer :: ptr_xVec(:), ptr_bVec(:), ptr_ResVec(:)
   PetscInt Lid,Gid,num_its, gqid,gfid,aid
   PetscInt ss,qq,ff
   PetscReal val, res_norm2, src_norm2, RHS_leading_coeff, Kss_coeff,inv_du(3)
   PetscErrorCode ierr
   integer istatus(MPI_STATUS_SIZE)
   PetscReal source(SG%dim_gq,SG%dim_a)
   PetscReal, allocatable :: xVec_mype(:),xVec_partd(:),bVec_mype(:),bVec_partd(:), ResVec_mype(:),ResVec_partd(:)
   real(wp) tt0,tt1
   real(wp) AD ! Adiabatic Factor = e*phi/Te0
   integer :: i,j,ij

   RHS_leading_coeff=1.602176634e-19 ! element_chargey (e)

   ! Assign and Assemble bVec (source term)
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      qq=mapG2P(Gid)%gqid
      ss=mapG2P(Gid)%aid
      ! val=RHS_leading_coeff*1e6*(SG%dni(qq,ss,1)-SG%dne(qq,ss,1))
         val=RHS_leading_coeff*1e6*(SG%dni_SM(qq,ss)-SG%dne_SM(qq,ss))
      call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   call Apply_BC_to_GK_Poisson_RHS_Coeffs


   call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
   call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)

   ! Solve the linear system
   tt0=MPI_WTIME()
   call KSPSolve(ksp_P,bVec,xVec,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0 .or. mype==1) write(stdout,*) "mype=",mype,"; KSPSolve takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call KSPView(ksp_P,PETSC_VIEWER_STDOUT_WORLD, ierr);CHKERRQ(ierr);

   ! Check the residual error and number of the iterations
   call MatMult(PMat,xVec,ResVec,ierr); CHKERRQ(ierr)
   call VecAXPY(ResVec,-1.0,bVec,ierr); CHKERRQ(ierr)
   call VecNorm(ResVec,NORM_2,res_norm2,ierr); CHKERRQ(ierr)
   call VecNorm(bVec,NORM_2,src_norm2,ierr); CHKERRQ(ierr)
   call KSPGetIterationNumber(ksp_P,num_its,ierr); CHKERRQ(ierr)
   if(mype==0) write(stdout,"(A,ES16.8,A,4I)") "Norm of error=", res_norm2/src_norm2, ", iterations=",num_its
   call flush(stdout)

   ! Copy xVec (solution) to the SG%test_phi vector
   allocate(xVec_mype(Lsize), xVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(xVec,ptr_xVec,ierr)
   xVec_mype(:)=ptr_xVec(1:Lsize)
   call VecRestoreArrayReadF90(xVec,ptr_xVec,ierr)
   call MPI_ALLGATHERV(xVec_mype,Lsize,mpi_Rsize,xVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%phi(mapG2P(Gid)%gqid,mapG2P(Gid)%aid,1) = xVec_partd(Lid)
   enddo
   SG%phi(SG%dim_gq,:,1)=SG%phi(1,:,1)
   deallocate(xVec_mype, xVec_partd)

   ! SG%phi(:,1:35,:) = 0.0
   ! do i=35,60
   !    SG%phi(:,i,:) = (i-35)*(1e4)/25.0
   ! enddo
   ! SG%phi(:,61:SG%dim_s-11,:) = 1e4
   ! do i=SG%dim_s-10,SG%dim_s
   !    SG%phi(:,i,:) = 1e4 -(i-(SG%dim_s-10))*(1e4)/10.0
   ! enddo

   ! do j=1,SG%dim_gq
   !    SG%phi(j,:,:) = SG%phi(j,:,:)*(1.0+0.5*cos(3*(j-1)/real(SG%dim_gq-1)*2.0*pi))
   ! enddo


   ! Copy bVec (source) to the SG%test_src vector
   allocate(bVec_mype(Lsize), bVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(bVec,ptr_bVec,ierr)
   bVec_mype(:)=ptr_bVec(1:Lsize)
   call VecRestoreArrayReadF90(bVec,ptr_bVec,ierr)
   call MPI_ALLGATHERV(bVec_mype,Lsize,mpi_Rsize,bVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%test_src_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = bVec_partd(Lid)
   enddo
   SG%test_src_2D(SG%dim_gq,:)=SG%test_src_2D(1,:)
   deallocate(bVec_mype, bVec_partd)

   ! Copy ResVec (residual=A*x-b) to the SG%test_res vector
   allocate(ResVec_mype(Lsize), ResVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(ResVec,ptr_ResVec,ierr)
   ResVec_mype(:)=ptr_ResVec(1:Lsize)
   call VecRestoreArrayReadF90(ResVec,ptr_ResVec,ierr)
   call MPI_ALLGATHERV(ResVec_mype,Lsize,mpi_Rsize,ResVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%test_res_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = ResVec_partd(Lid)
   enddo
   SG%test_res_2D(SG%dim_gq,:)=SG%test_res_2D(1,:)
   deallocate(ResVec_mype, ResVec_partd)

   call MPI_SENDRECV(SG%phi(:,:,1),SG%dim_gq*SG%dim_a,mpi_Rsize,right_pe,0,&
                     SG%phi(:,:,0),SG%dim_gq*SG%dim_a,mpi_Rsize,left_pe,0,toroidal_comm,istatus,ierr)
   call MPI_SENDRECV(SG%phi(:,:,1),SG%dim_gq*SG%dim_a,mpi_Rsize,left_pe,0,&
                     SG%phi(:,:,2),SG%dim_gq*SG%dim_a,mpi_Rsize,right_pe,0,toroidal_comm,istatus,ierr)


   ! Copy SG%phi to SG_phi that will be used on GPU side
   SG_phi = SG%phi
   !$acc update device(SG_phi)

   if(mod(istep,msnap)==0) then
      ! To get full 3D data for diagnostic purpose
      call MPI_ALLGATHER(SG%phi(:,:,1),SG%dim_gq*SG%dim_a,mpi_Rsize,SG%test_phi_3D,SG%dim_gq*SG%dim_a,mpi_Rsize,toroidal_comm,ierr)
      call MPI_ALLGATHER(SG%test_src_2D,SG%dim_gq*SG%dim_a,mpi_Rsize,SG%test_src_3D,SG%dim_gq*SG%dim_a,mpi_Rsize,toroidal_comm,ierr)
   endif

end subroutine Solve_Poisson3D_on_SG

subroutine Solve_Nonlinear_Poisson3D_on_SG
   use precision, only:wp
   ! use global_parameters, only:stdout,mype,ifelectron,mzetamax,mgrid
   use global_parameters
   use field_array

   use Structured_Grid
   use particle_decomp
   implicit none

   PetscScalar, pointer :: ptr_xVec(:), ptr_srcVec(:), ptr_ResVec(:)
   PetscInt Lid,Gid,num_its, gqid,gfid,aid
   PetscInt ss,qq,ff
   PetscReal val,src0, res_norm2, src_norm2, ELEMENT_CHARGE
   PetscReal delta_Phi_norm2, ori_Phi_norm2, RelErr
   PetscErrorCode ierr
   integer istatus(MPI_STATUS_SIZE)
   PetscReal source(SG%dim_gq,SG%dim_a)
   PetscReal, allocatable :: xVec_mype(:),xVec_partd(:),srcVec_mype(:),srcVec_partd(:), ResVec_mype(:),ResVec_partd(:)
   real(wp) tt0,tt1
   real(wp) AD ! Adiabatic Factor = e*phi/Te0
   integer :: i,j,ij,loop

   ELEMENT_CHARGE=1.602176634e-19 ! element_chargey (e)


   do loop=1,10
      !ptr_xVec is pointing xVec which is Phi
      call VecGetArrayReadF90(xVec,ptr_xVec,ierr) !xVec is Phi
      ! Assign and Assemble bVec (source term)
      do Lid=1,Lsize
         Gid=mapL2G(Lid)
         qq=mapG2P(Gid)%gqid
         ss=mapG2P(Gid)%aid

         ! RHS of the Poisson Equation (src term): e(ni0+dni-dne)
         src0 = ELEMENT_CHARGE*1e6*(SG%n0_1D(ss)+SG%dni_SM(qq,ss)-SG%dne_SM(qq,ss))
         call VecSetValues(srcVec,1,Gid,src0,INSERT_VALUES,ierr); CHKERRQ(ierr)

         ! Second term of Fvec: e*ne0*exp(e*phi/Te0)
         val = ELEMENT_CHARGE*1e6*SG%n0_1D(ss)*exp(ptr_xVec(Lid)/SG%Te0_eV_1D(ss))

         ! Add third term (which is -src0)
         val = val -src0
         call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      enddo
      call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
      call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)

      ! Calculate bVec = PMat*xVec+bVec
      ! So, bVec = A*phi-e(n0*(1-exp(AD)+dni-dne)
      call MatMultAdd(PMat,xVec,bVec,bVec,ierr); CHKERRQ(ierr)
      call Apply_BC_to_GK_Poisson_RHS_Coeffs
      call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
      call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)
      call VecAssemblyBegin(srcVec,ierr); CHKERRQ(ierr)
      call VecAssemblyEnd(srcVec,ierr); CHKERRQ(ierr)


      ! Update Jmat
      ! Set PMat value as the first term
      call MatCopy(PMat,Jmat,SAME_NONZERO_PATTERN,ierr); CHKERRQ(ierr)
      ! Add another term from the nonlinear adiabatic response
      do Lid=1,Lsize
         Gid=mapL2G(Lid)
         ss = mapG2P(Gid)%aid
         val = 1e6*SG%n0_1D(ss)*ELEMENT_CHARGE/SG%Te0_eV_1D(ss)*exp(ptr_xVec(Lid)/SG%Te0_eV_1D(ss))
         call MatSetValues(Jmat,1,Gid,1,Gid,val,ADD_VALUES,ierr); CHKERRQ(ierr)
      enddo
      ! Assembly Jacobian Matrix
      tt0=MPI_WTIME()
      call MatAssemblyBegin(Jmat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(Jmat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "JMat Assembly takes ",tt1-tt0," secs"
      call flush(stdout)

      ! Update pre_Jmat
      if(loop==1 .and. irk==1) then
         call MatCopy(Jmat,pre_Jmat,SAME_NONZERO_PATTERN,ierr); CHKERRQ(ierr)
         if(mype==0) write(stdout,*) "Updated pre_Jmat!"
         call flush(stdout)
      endif


      ! Restore pointer
      call VecRestoreArrayReadF90(xVec,ptr_xVec,ierr)


      ! Solve the linear system (J*y=b; so y=J^(-1)*b)
      tt0=MPI_WTIME()
      call KSPSolve(ksp_P,bVec,yVec,ierr); CHKERRQ(ierr)
      tt1=MPI_WTIME()
      if(mype==0 .or. mype==1) write(stdout,*) "mype=",mype,"; KSPSolve takes ",tt1-tt0," secs"
      call flush(stdout)

      ! Check the residual error and number of the iterations
      call MatMult(Jmat,yVec,ResVec,ierr); CHKERRQ(ierr)
      call VecAXPY(ResVec,-1.0,bVec,ierr); CHKERRQ(ierr)
      call VecNorm(ResVec,NORM_2,res_norm2,ierr); CHKERRQ(ierr)
      call VecNorm(bVec,NORM_2,src_norm2,ierr); CHKERRQ(ierr)
      call KSPGetIterationNumber(ksp_P,num_its,ierr); CHKERRQ(ierr)
      if(mype==0) write(stdout,"(A,ES16.8,A,4I)") "Solving y=J^(-1)*F; Norm of error=", res_norm2/src_norm2, ", iterations=",num_its
      call flush(stdout)


      call VecNorm(yVec,NORM_2,delta_Phi_norm2,ierr); CHKERRQ(ierr)
      call VecNorm(xVec,NORM_2,ori_Phi_norm2,ierr); CHKERRQ(ierr)
      ! Calculate new potential of (k+1)-th step: xVec(k+1)=xVec(k)-yVec
      call VecAXPY(xVec,-1.0,yVec,ierr); CHKERRQ(ierr)


      ! Calculate |Residual|/|Src| = norm2(bVec)/norm2(srcVec)
      call VecNorm(bVec,NORM_2,res_norm2,ierr); CHKERRQ(ierr)
      call VecNorm(srcVec,NORM_2,src_norm2,ierr); CHKERRQ(ierr)
      RelErr=res_norm2/src_norm2
      if(mype==0) then
         write(stdout,"(A,I4,A,ES16.8)") "loop=",loop, "; Relative Error(|Res|/|Src|)=", RelErr
         write(stdout,"(A,I4,A,ES16.8)") "loop=",loop, "; evolving phi ratio=", delta_Phi_norm2/ori_Phi_norm2
         call flush(stdout)
      endif

      if(RelErr<1e-9 .or. delta_Phi_norm2/ori_Phi_norm2<1e-7) then
         if(mype==0) then
            write(stdout,"(A,I4,A,ES16.8)") "last_loop=",loop, "; Relative Error(|Res|/|Src|)=", RelErr
            write(stdout,"(A,I4,A,ES16.8)") "last_loop=",loop, "; evolving phi ratio=", delta_Phi_norm2/ori_Phi_norm2
            call flush(stdout)
         endif
         exit
      endif
   enddo

   ! Copy xVec (solution) to the SG%test_phi vector
   allocate(xVec_mype(Lsize), xVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(xVec,ptr_xVec,ierr)
   xVec_mype(:)=ptr_xVec(1:Lsize)
   call VecRestoreArrayReadF90(xVec,ptr_xVec,ierr)
   call MPI_ALLGATHERV(xVec_mype,Lsize,mpi_Rsize,xVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%phi(mapG2P(Gid)%gqid,mapG2P(Gid)%aid,1) = xVec_partd(Lid)
   enddo
   SG%phi(SG%dim_gq,:,1)=SG%phi(1,:,1)
   deallocate(xVec_mype, xVec_partd)

   ! Communicate phi in toroidal direction
   call MPI_SENDRECV(SG%phi(:,:,1),SG%dim_gq*SG%dim_a,mpi_Rsize,right_pe,0,&
                     SG%phi(:,:,0),SG%dim_gq*SG%dim_a,mpi_Rsize,left_pe,0,toroidal_comm,istatus,ierr)
   call MPI_SENDRECV(SG%phi(:,:,1),SG%dim_gq*SG%dim_a,mpi_Rsize,left_pe,0,&
                     SG%phi(:,:,2),SG%dim_gq*SG%dim_a,mpi_Rsize,right_pe,0,toroidal_comm,istatus,ierr)


   ! Copy SG%phi to SG_phi that will be used on GPU side
   SG_phi = SG%phi
   !$acc update device(SG_phi)


   if(mod(istep,msnap)==0) then
      ! Copy srcVec (source) to the SG%test_src vector
      allocate(srcVec_mype(Lsize), srcVec_partd(sum_Lsize_partd))
      call VecGetArrayReadF90(srcVec,ptr_srcVec,ierr)
      srcVec_mype(:)=ptr_srcVec(1:Lsize)
      call VecRestoreArrayReadF90(srcVec,ptr_srcVec,ierr)
      call MPI_ALLGATHERV(srcVec_mype,Lsize,mpi_Rsize,srcVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
      do Lid=1,sum_Lsize_partd
         Gid=recv_mapL2G_partd(Lid)
         SG%test_src_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = srcVec_partd(Lid)
      enddo
      SG%test_src_2D(SG%dim_gq,:)=SG%test_src_2D(1,:)
      deallocate(srcVec_mype, srcVec_partd)


      call VecCopy(bVec,ResVec,ierr); CHKERRQ(ierr)
      ! Copy ResVec (residual=A*x-b) to the SG%test_res vector
      allocate(ResVec_mype(Lsize), ResVec_partd(sum_Lsize_partd))
      call VecGetArrayReadF90(ResVec,ptr_ResVec,ierr)
      ResVec_mype(:)=ptr_ResVec(1:Lsize)
      call VecRestoreArrayReadF90(ResVec,ptr_ResVec,ierr)
      call MPI_ALLGATHERV(ResVec_mype,Lsize,mpi_Rsize,ResVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
      do Lid=1,sum_Lsize_partd
         Gid=recv_mapL2G_partd(Lid)
         SG%test_res_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = ResVec_partd(Lid)
      enddo
      SG%test_res_2D(SG%dim_gq,:)=SG%test_res_2D(1,:)
      deallocate(ResVec_mype, ResVec_partd)


   endif

end subroutine Solve_Nonlinear_Poisson3D_on_SG


subroutine Solve_test_Poisson3D_on_SG()
   use precision, only:wp
   ! use global_parameters, only:stdout,mype,ifelectron,mzetamax,mgrid
   use global_parameters
   use field_array

   use Structured_Grid
   use particle_decomp
   use field_array
   implicit none
   PetscScalar, pointer :: ptr_xVec(:), ptr_bVec(:), ptr_ResVec(:)
   PetscInt Lid,Gid,num_its, gqid,gfid,aid
   PetscInt ss,qq,ff
   PetscReal val, res_norm2, src_norm2, RHS_leading_coeff, Kss_coeff,inv_du(3)
   PetscErrorCode ierr
   PetscReal source(SG%dim_gq,SG%dim_a)
   PetscReal, allocatable :: xVec_mype(:),xVec_partd(:),bVec_mype(:),bVec_partd(:), ResVec_mype(:),ResVec_partd(:)
   real(wp) tt0,tt1
   real(wp), external :: density,Te0_eV
   real(wp) source_from_GTS(mgrid), dne_adiab_ratio, intpl_dphi_NZ
   integer :: i,j,ij

   ! source = SG%iScaData(:,:,1,1)-(SG%eScaData(:,:,1,1)+int(ifelectron/2)*SG%dne_adiab)
   source=0.0
   ! do aid=20,30
   !    do gqid=1,SG%dim_gq
   !       source(gqid,aid)=sin(6*(SG%gq_1D(gqid))-4*(SG%zetamin))
   !    enddo
   ! enddo
   ! call Interpolate_OriGridData_into_SG(test_src,source)

   ! do i=0,mpsi
   !    do j=0,mtheta(i)
   !       ij=igrid(i) + j

   !       intpl_dphi_NZ = phi(1,ij) -real(1.0-ZF_SEPARATION)*phi00(i)
   !       dne_adiab_ratio = abs(csc0)*intpl_dphi_NZ/Te_grid(i)
   !       source_from_GTS(ij)=(densityi(1,ij)+real(ZF_SEPARATION)*zonali(i)) - (densitye(1,ij)+real(ZF_SEPARATION)*zonale(i) + int(ifelectron/2)*dne_adiab_ratio)
   !       source_from_GTS(ij)=source_from_GTS(ij) + int(ifelectron/2)*abs(csc0)*(phi(1,ij)+real(ZF_SEPARATION)*phi00(i))/Te_grid(i)
   !    enddo
   ! enddo
   ! call Interpolate_OriGridData_into_SG(source_from_GTS,source)


   source=SG%test_src_2D
   ! do aid=20,25
   !    do gqid=1,SG%dim_gq
   !       source(gqid,aid)=0.02
   !    enddo
   ! enddo


   RHS_leading_coeff=1.602176634e-19 ! element_charge/permittivity (e/eps0)

   source=RHS_leading_coeff*source
   do aid=1,SG%dim_a
      source(:,aid) = 1e6*density(SG%a_1D(aid),0.0)*source(:,aid)
   enddo

   !add the adiabatic electron response term in the LHS (PMat)
   ! NOTE: The following line must locate after "call Apply_BC_to_GK_Poisson_LHS_Coeffs"
   if(ifelectron==1 .and. Flag_Poisson3D_on_SG==2) then
      ! To mimic GTS Way
      if(adiab > 0.0_wp) then
         do Lid=1,Lsize
            Gid= mapL2G(Lid)
            qq = mapG2P(Gid)%gqid
            ss = mapG2P(Gid)%aid
            source(qq,ss) = source(qq,ss) + adiab*(1.602176634e-19)*1e6*density(SG%a_1D(ss),SG%gq_1D(qq))/Te0_eV(SG%a_1D(ss),SG%gq_1D(qq))*SG%test_phi_2D(qq,ss)
         enddo
      else
         do Lid=1,Lsize
            Gid= mapL2G(Lid)
            qq = mapG2P(Gid)%gqid
            ss = mapG2P(Gid)%aid
            source(qq,ss) = source(qq,ss) + 0.01/(Cs*Cs)/(1e3*Tem_00)*(1.602176634e-19)*1e6*density(SG%a_1D(ss),SG%gq_1D(qq))*SG%test_phi_2D(qq,ss)
         enddo
      endif
   endif

   ! Assign and Assemble bVec (source term)
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      val=source(mapG2P(Gid)%gqid,mapG2P(Gid)%aid)
      call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   call Apply_BC_to_GK_Poisson_RHS_Coeffs




   call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
   call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)

   ! call VecView(bVec,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   ! Solve the linear system
   tt0=MPI_WTIME()
   call KSPSolve(ksp_P,bVec,xVec,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPSolve takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call KSPView(ksp_P,PETSC_VIEWER_STDOUT_WORLD, ierr);CHKERRQ(ierr);

   ! Check the residual error and number of the iterations
   call MatMult(PMat,xVec,ResVec,ierr); CHKERRQ(ierr)
   call VecAXPY(ResVec,-1.0,bVec,ierr); CHKERRQ(ierr)
   call VecNorm(ResVec,NORM_2,res_norm2,ierr); CHKERRQ(ierr)
   call VecNorm(bVec,NORM_2,src_norm2,ierr); CHKERRQ(ierr)
   call KSPGetIterationNumber(ksp_P,num_its,ierr); CHKERRQ(ierr)
   if(mype==0) write(stdout,"(A,ES16.8,A,4I)") "Norm of error=", res_norm2/src_norm2, ", iterations=",num_its
   call flush(stdout)

   ! Copy xVec (solution) to the SG%test_phi vector
   allocate(xVec_mype(Lsize), xVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(xVec,ptr_xVec,ierr)
   xVec_mype(:)=ptr_xVec(1:Lsize)
   call VecRestoreArrayReadF90(xVec,ptr_xVec,ierr)
   call MPI_ALLGATHERV(xVec_mype,Lsize,mpi_Rsize,xVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%test_phi_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = xVec_partd(Lid)
   enddo
   SG%test_phi_2D(SG%dim_gq,:)=SG%test_phi_2D(1,:)
   deallocate(xVec_mype, xVec_partd)


   ! Copy bVec (source) to the SG%test_phi vector
   allocate(bVec_mype(Lsize), bVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(bVec,ptr_bVec,ierr)
   bVec_mype(:)=ptr_bVec(1:Lsize)
   call VecRestoreArrayReadF90(bVec,ptr_bVec,ierr)
   call MPI_ALLGATHERV(bVec_mype,Lsize,mpi_Rsize,bVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%test_src_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = bVec_partd(Lid)
   enddo
   SG%test_src_2D(SG%dim_gq,:)=SG%test_src_2D(1,:)
   deallocate(bVec_mype, bVec_partd)

   ! Copy ResVec (residual=A*x-b) to the SG%test_phi vector
   allocate(ResVec_mype(Lsize), ResVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(ResVec,ptr_ResVec,ierr)
   ResVec_mype(:)=ptr_ResVec(1:Lsize)
   call VecRestoreArrayReadF90(ResVec,ptr_ResVec,ierr)
   call MPI_ALLGATHERV(ResVec_mype,Lsize,mpi_Rsize,ResVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%test_res_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = ResVec_partd(Lid)
   enddo
   SG%test_res_2D(SG%dim_gq,:)=SG%test_res_2D(1,:)
   deallocate(ResVec_mype, ResVec_partd)



end subroutine Solve_test_Poisson3D_on_SG


subroutine Cal_IonPolarDensity_and_RealCharge
   use precision, only:wp
   use global_parameters
   use field_array
   use Structured_Grid
   use particle_decomp
   implicit none

   PetscScalar, pointer :: ptr_Vec(:)
   PetscReal, allocatable :: Vec_mype(:), Vec_partd(:)
   integer Lid,Gid
   Vec RHS_Vec ! Right Hand Side of PMat*xVec
   PetscErrorCode ierr

   call VecCreate(PETSC_COMM_WORLD,RHS_Vec,ierr); CHKERRQ(ierr)
   call VecSetSizes(RHS_Vec,Lsize,Gsize,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(RHS_Vec,ierr); CHKERRQ(ierr)

   allocate(Vec_mype(Lsize), Vec_partd(sum_Lsize_partd))


   ! Calculate LHS of the pure Poission Equation (not GK one) to get the real charge density
   ! Note that xVec is the potential calculated already from other routine
   ! Here, RHS vector is the Real Charge Density (including Ion Polarization Density) = ni_GC+ni_Pol-ne
   call Update_GK_Poisson_LHS_Coeffs(1) ! For Pure Poisson Equation
   call MatMult(PMat,xVec,RHS_Vec,ierr); CHKERRQ(ierr)
   call VecGetArrayReadF90(RHS_Vec,ptr_Vec,ierr)
   Vec_mype(:)=ptr_Vec(1:Lsize)
   call VecRestoreArrayReadF90(RHS_Vec,ptr_Vec,ierr)
   call MPI_ALLGATHERV(Vec_mype,Lsize,mpi_Rsize,Vec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%RealCharge_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = Vec_partd(Lid)
   enddo
   SG%RealCharge_2D(SG%dim_gq,:)=SG%RealCharge_2D(1,:)


   ! Calculate Ion Polarization Density operator
   ! Note the signs of the equation
   ! niPol=(nabla*Pol_T*nabla*phi)/e=(-PMat*xVec)/e
   ! Here, RHS is the Ion Polarization Density
   call Update_GK_Poisson_LHS_Coeffs(2) ! For Ion Polarization Density
   call MatMult(PMat,xVec,RHS_Vec,ierr); CHKERRQ(ierr)
   call VecGetArrayReadF90(RHS_Vec,ptr_Vec,ierr)
   Vec_mype(:)=ptr_Vec(1:Lsize)
   call VecRestoreArrayReadF90(RHS_Vec,ptr_Vec,ierr)
   call MPI_ALLGATHERV(Vec_mype,Lsize,mpi_Rsize,Vec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      SG%niPol_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = Vec_partd(Lid)
   enddo
   SG%niPol_2D(SG%dim_gq,:)=SG%niPol_2D(1,:)
   SG%niPol_2D=SG%niPol_2D/(-1.602176634e-19) ! Convert from charge denstiy to particle density

   ! Free the pointer and array
   call VecDestroy(RHS_Vec,ierr); CHKERRQ(ierr)
   deallocate(Vec_mype,Vec_partd)

   if(mod(istep,msnap)==0) then
      ! To get full 3D data for diagnostic purpose
      call MPI_ALLGATHER(SG%RealCharge_2D,SG%dim_gq*SG%dim_a,mpi_Rsize,SG%RealCharge_3D,SG%dim_gq*SG%dim_a,mpi_Rsize,toroidal_comm,ierr)
      call MPI_ALLGATHER(SG%niPol_2D,SG%dim_gq*SG%dim_a,mpi_Rsize,SG%niPol_3D,SG%dim_gq*SG%dim_a,mpi_Rsize,toroidal_comm,ierr)
   endif


   ! Recover original PMat (which is for GK Poisson Equation)
   call Update_GK_Poisson_LHS_Coeffs(0)

end subroutine Cal_IonPolarDensity_and_RealCharge

subroutine test_Laplacian_2D()
#include <petsc/finclude/petsc.h>
   use petscksp
   use global_parameters,only:stdout,mype,numberpe
   implicit none

   Vec x,b,u
   Mat PMat
   KSP ksp_P
   PetscReal val,x_norm2
   PetscInt :: gid, fGidx, MGidx_start, MGidx_end, VGidx_start, VGidx_end, nx=30,ny=30, neq, idxm, num_its
   PetscErrorCode ierr

   real tt0,tt1

#if defined(PETSC_USE_LOG)
   PetscLogStage stage;
#endif


   neq=nx*ny


   ! call PetscInitialize('./PETSc_input',ierr)
   ! if(ierr) return


   if(mype==0) write(stdout,*) "numberpe=",numberpe
   call flush(stdout)

   call MatCreate(PETSC_COMM_WORLD,PMat,ierr); CHKERRQ(ierr)
   call MatSetSizes(PMat,PETSC_DECIDE,PETSC_DECIDE,neq,neq,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(PMat,ierr); CHKERRQ(ierr)
   ! Note that the following d_nz=5, o_nz=5 are not optimized values.
   call MatMPIAIJSetPreallocation(PMat,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,ierr); CHKERRQ(ierr)
   call MatGetOwnershipRange(PMat,MGidx_start,MGidx_end,ierr); CHKERRQ(ierr)

   MGidx_end=MGidx_end-1 !Note that this is an global index of the last local row more proper for Fortran do loop

   call PetscLogStageRegister("Assembly", stage,ierr); CHKERRQ(ierr)
   call PetscLogStagePush(stage,ierr); CHKERRQ(ierr)


      ! do gid=0,neq-1
      !    idxm=gid-1; val=-1.0
      !    if(idxm>=0) call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      !    idxm=gid+1; val=-1.0
      !    if(idxm<=neq-1) call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      !    idxm=gid-ny; val=-1.0
      !    if(idxm>=0) call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      !    idxm=gid+ny; val=-1.0
      !    if(idxm<=neq-1) call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)

      !    idxm=gid; val=4.0
      !    call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr)
      ! enddo

   do gid=MGidx_start,MGidx_end
      idxm=gid-1; val=-1.0
      if(idxm>=0) call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      idxm=gid+1; val=-1.0
      if(idxm<=neq-1) call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      idxm=gid-ny; val=-1.0
      if(idxm>=0) call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      idxm=gid+ny; val=-1.0
      if(idxm<=neq-1) call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)

      idxm=gid; val=4.0
      call MatSetValues(PMat,1,gid,1,idxm,val,INSERT_VALUES,ierr)
   enddo
   tt0=MPI_WTIME()
   call MatAssemblyBegin(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)
   call PetscLogStagePop(ierr); CHKERRQ(ierr)

   ! call MatView(PMat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

   tt0=MPI_WTIME()
   call VecCreate(PETSC_COMM_WORLD,u,ierr); CHKERRQ(ierr)
   call VecSetSizes(u,PETSC_DECIDE,neq,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(u,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "VecCreate u takes ",tt1-tt0," secs"
   call flush(stdout)

   call VecDuplicate(u,b,ierr); CHKERRQ(ierr)
   call VecDuplicate(b,x,ierr); CHKERRQ(ierr)

   call VecGetOwnershipRange(u,VGidx_start,VGidx_end,ierr); CHKERRQ(ierr)
   VGidx_end=VGidx_end-1 !Note that this is an global index of the last local row more proper for Fortran do loop
   ! write(stdout,"(5(A,1x,I8))")"mype=",mype,", MGidx_start=",MGidx_start,", MGidx_end=",MGidx_end,", VGidx_start=",VGidx_start, ", VGidx_end=",VGidx_end

   call VecSet(u,1.0,ierr); CHKERRQ(ierr)

   ! call VecView(u,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   tt0=MPI_WTIME()
   call MatMult(PMat,u,b,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatMult takes ",tt1-tt0," secs"
   call flush(stdout)


   tt0=MPI_WTIME()
   call KSPCreate(PETSC_COMM_WORLD,ksp_P,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp_P,PMat,PMat,ierr); CHKERRQ(ierr)
   ! call KSPSetTolerances(ksp_P,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp_P,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
   call flush(stdout)

   tt0=MPI_WTIME()
   ! Solve the linear system
   call KSPSolve(ksp_P,b,x,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPSolve takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   ! Check solution and error
   call VecAXPY(x,-1.0,u,ierr); CHKERRQ(ierr)
   call VecNorm(x,NORM_2,x_norm2,ierr); CHKERRQ(ierr)
   call KSPGetIterationNumber(ksp_P,num_its,ierr); CHKERRQ(ierr)

   if(mype==0) write(stdout,"(A,ES16.8,A,4I)") "Norm of error=", x_norm2, ", iterations=",num_its

   call MPI_Barrier(PETSC_COMM_WORLD,ierr)
   write(stdout,"(A,I8,A,(I8,A,I8),A,(I8,A,I8),A)")"mype=",mype,", idx_start_(M,V)=(",MGidx_start,",",VGidx_start,"), idx_end_(M,V)=(",MGidx_end,",",VGidx_end,")"
   call flush(stdout)

   call KSPDestroy(ksp_P,ierr); CHKERRQ(ierr)
   call VecDestroy(u,ierr); CHKERRQ(ierr)
   call VecDestroy(x,ierr); CHKERRQ(ierr)
   call VecDestroy(b,ierr); CHKERRQ(ierr)
   call MatDestroy(PMat,ierr); CHKERRQ(ierr)

   ! call PetscFinalize(ierr)

   return
end subroutine test_Laplacian_2D


subroutine Construct_GK_Poisson3D_Solver_using_Structured_Grid

   use Structured_Grid
   use particle_decomp
   use global_parameters,only:mype,stdout,zetamax,zetamin,mzetamax
   implicit none

   PetscInt Lid,Gid,qq,ss,ff,tnum
   PetscReal off_val
   real tt0, tt1
   PetscReal inv_Jacob
   PetscErrorCode ierr
   real, external:: dsda
   integer :: dim_gq, dim_a


   dim_a =SG%dim_a
   dim_gq=SG%dim_gq

   ! call Initialize_MGY_PETSc_Interface_Data(dim_gq,dim_a,ntoroidal)

   ! Tensors on the poloidal plane, and Coeff Tensor including neighboring left and right poloidal plane
   ! Note that there are ghost cells for easier indexing
   allocate(Polar_T(3,3,0:dim_gq,0:dim_a+1))
   allocate(Eps_T(3,3,0:dim_gq,0:dim_a+1,-1:1), CT(3,3,0:dim_gq,0:dim_a+1,-1:1))
   call Determine_GK_PoissonTensors_on_Structured_Grid(0)

   max_nnz=19 ! maximum number of non-zero values per each row
   allocate(tarGid(max_nnz,Lsize), tarVal(max_nnz,Lsize))

   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   do Lid=1,Lsize
      Gid=mapL2G(Lid)

      qq = mapG2P(Gid)%gqid
      ss = mapG2P(Gid)%aid
      ff = mapG2P(Gid)%gfid

      inv_Jacob=1.0/SG%Jacob_sqf(qq,ss)


      ! 2 (gqid-1) & 3 (gqid+1)
      off_val        =0.125*(CT(1,2,qq,ss+1,0)-CT(1,2,qq,ss-1,0) +CT(3,2,qq,ss,+1)-CT(3,2,qq,ss,-1))
      tarGid(2 ,Lid) =mapP2G(qq-1,ss,ff)
      tarVal(2 ,Lid) =0.5*(CT(2,2,qq,ss,0)+CT(2,2,qq-1,ss,0)) - off_val
      tarGid(3 ,Lid) =mapP2G(qq+1,ss,ff)
      tarVal(3 ,Lid) =0.5*(CT(2,2,qq,ss,0)+CT(2,2,qq+1,ss,0)) + off_val

      ! 4 (sid-1) & 5 (sid+1)
      off_val        =0.125*(CT(2,1,qq+1,ss,0)-CT(2,1,qq-1,ss,0) +CT(3,1,qq,ss,+1)-CT(3,1,qq,ss,-1))
      tarGid(4 ,Lid) =mapP2G(qq,ss-1,ff)
      tarVal(4 ,Lid) =0.5*(CT(1,1,qq,ss,0)+CT(1,1,qq,ss-1,0)) - off_val
      tarGid(5 ,Lid) =mapP2G(qq,ss+1,ff)
      tarVal(5 ,Lid) =0.5*(CT(1,1,qq,ss,0)+CT(1,1,qq,ss+1,0)) + off_val

      ! 6 (gfid-1) & 7 (gfid+1)
      off_val        =0.125*(CT(1,3,qq,ss+1,0)-CT(1,3,qq,ss-1,0) +CT(2,3,qq+1,ss,0)-CT(2,3,qq-1,ss,0))
      tarGid(6 ,Lid) =mapP2G(qq,ss,ff-1)
      tarVal(6 ,Lid) =0.5*(CT(3,3,qq,ss,0)+CT(3,3,qq,ss,-1)) - off_val
      tarGid(7 ,Lid) =mapP2G(qq,ss,ff+1)
      tarVal(7 ,Lid) =0.5*(CT(3,3,qq,ss,0)+CT(3,3,qq,ss,+1)) + off_val


      ! 8 (sid+1,qid+1) & 9 (sid-1,qid-1)
      tarGid(8 ,Lid) =mapP2G(qq+1,ss+1,ff)
      tarVal(8 ,Lid) =+0.125*(CT(1,2,qq,ss,0)+CT(1,2,qq,ss+1,0)+CT(2,1,qq,ss,0)+CT(2,1,qq+1,ss,0))
      tarGid(9 ,Lid) =mapP2G(qq-1,ss-1,ff)
      tarVal(9 ,Lid) =+0.125*(CT(1,2,qq,ss,0)+CT(1,2,qq,ss-1,0)+CT(2,1,qq,ss,0)+CT(2,1,qq-1,ss,0))
      ! 10 (sid+1,qid-1) & 11 (sid-1,qid+1)
      tarGid(10,Lid) =mapP2G(qq-1,ss+1,ff)
      tarVal(10,Lid) =-0.125*(CT(1,2,qq,ss,0)+CT(1,2,qq,ss+1,0)+CT(2,1,qq,ss,0)+CT(2,1,qq-1,ss,0))
      tarGid(11,Lid) =mapP2G(qq+1,ss-1,ff)
      tarVal(11,Lid) =-0.125*(CT(1,2,qq,ss,0)+CT(1,2,qq,ss-1,0)+CT(2,1,qq,ss,0)+CT(2,1,qq+1,ss,0))

      ! 12 (sid+1,gfid+1) & 13 (sid-1,gfid-1)
      tarGid(12,Lid) =mapP2G(qq,ss+1,ff+1)
      tarVal(12,Lid) =+0.125*(CT(1,3,qq,ss,0)+CT(1,3,qq,ss+1,0)+CT(3,1,qq,ss,0)+CT(3,1,qq,ss,+1))
      tarGid(13,Lid) =mapP2G(qq,ss-1,ff-1)
      tarVal(13,Lid) =+0.125*(CT(1,3,qq,ss,0)+CT(1,3,qq,ss-1,0)+CT(3,1,qq,ss,0)+CT(3,1,qq,ss,-1))
      ! 14 (sid+1,gfid-1) & 15 (sid-1,gfid+1)
      tarGid(14,Lid) =mapP2G(qq,ss+1,ff-1)
      tarVal(14,Lid) =-0.125*(CT(1,3,qq,ss,0)+CT(1,3,qq,ss+1,0)+CT(3,1,qq,ss,0)+CT(3,1,qq,ss,-1))
      tarGid(15,Lid) =mapP2G(qq,ss-1,ff+1)
      tarVal(15,Lid) =-0.125*(CT(1,3,qq,ss,0)+CT(1,3,qq,ss-1,0)+CT(3,1,qq,ss,0)+CT(3,1,qq,ss,+1))

      ! 16 (gqid+1,gfid+1) & 17  (gqid-1, gfid-1)
      tarGid(16,Lid) =mapP2G(qq+1,ss,ff+1)
      tarVal(16,Lid) =+0.125*(CT(2,3,qq,ss,0)+CT(2,3,qq+1,ss,0)+CT(3,2,qq,ss,0)+CT(3,2,qq,ss,+1))
      tarGid(17,Lid) =mapP2G(qq-1,ss,ff-1)
      tarVal(17,Lid) =+0.125*(CT(2,3,qq,ss,0)+CT(2,3,qq-1,ss,0)+CT(3,2,qq,ss,0)+CT(3,2,qq,ss,-1))
      ! 18 (gqid+1,gfid-1) & 18 (gqid-1, gfid+1)
      tarGid(18,Lid) =mapP2G(qq+1,ss,ff-1)
      tarVal(18,Lid) =-0.125*(CT(2,3,qq,ss,0)+CT(2,3,qq+1,ss,0)+CT(3,2,qq,ss,0)+CT(3,2,qq,ss,-1))
      tarGid(19,Lid) =mapP2G(qq-1,ss,ff+1)
      tarVal(19,Lid) =-0.125*(CT(2,3,qq,ss,0)+CT(2,3,qq-1,ss,0)+CT(3,2,qq,ss,0)+CT(3,2,qq,ss,+1))


      ! 1 (gqid,sid,gfid)
      tarGid(1,Lid)=mapP2G(qq,ss,ff)
      tarVal(1,Lid)=-sum(tarVal(2:19,Lid))

      ! IMPORTANT!! Multiply -inv_Jacob to the above coefficients,
      ! so make sure that target coefficient values represent (-nabla*(eps_T*nabla*phi))
      tarVal(:,Lid) = -1.0*inv_Jacob*tarVal(:,Lid)

   enddo

   call Apply_BC_to_GK_Poisson_LHS_Coeffs


   ! To calculate numbers of non-zero elements for diagonal(d_nnz) & off-diagonal(o_nnz) blocks
   allocate(d_nnz(Lsize), o_nnz(Lsize))
   d_nnz=0; o_nnz=0;
   do Lid=1,Lsize
      do tnum=1,max_nnz
         if(tarGid(tnum,Lid)>=0) then
            if(tarGid(tnum,Lid)<Gid_start .or. tarGid(tnum,Lid)>Gid_end) then
               o_nnz(Lid)=o_nnz(Lid)+1
            else
               d_nnz(Lid)=d_nnz(Lid)+1
            endif
         endif
      enddo
   enddo


   ! Construct Matrix (PMat)
   call MatCreate(PETSC_COMM_WORLD,PMat,ierr); CHKERRQ(ierr)
   call MatSetSizes(PMat,Lsize,Lsize,Gsize,Gsize,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(PMat,ierr); CHKERRQ(ierr)
   ! Preallocate the memory using o_nnz and d_nnz
   call MatMPIAIJSetPreallocation(PMat,max_nnz,d_nnz,max_nnz,o_nnz,ierr); CHKERRQ(ierr)
   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(PMat,1,Gid,max_nnz,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call MatView(PMat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)


   call VecCreate(PETSC_COMM_WORLD,bVec,ierr); CHKERRQ(ierr)
   call VecSetSizes(bVec,Lsize,Gsize,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(bVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,xVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,ResVec,ierr); CHKERRQ(ierr)

   ! Construct KSP
   tt0=MPI_WTIME()
   call KSPCreate(PETSC_COMM_WORLD,ksp_P,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp_P,PMat,PMat,ierr); CHKERRQ(ierr)
   ! call KSPSetTolerances(ksp_P,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp_P,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
   call flush(stdout)


end subroutine Construct_GK_Poisson3D_Solver_using_Structured_Grid

subroutine Construct_GK_Poisson3D_Solver_with_Extremal_Method

   use Structured_Grid
   use particle_decomp
   use global_parameters
   use field_array, only:FLA_Te0
   implicit none

   PetscInt Lid,Gid,qq,ss,ff,tnum
   real tt0, tt1
   PetscReal inv_Jacob
   PetscErrorCode ierr
   real, external:: dsda
   integer :: dim_gq, dim_a
   PetscReal AntiSym1, AntiSym2, Sym1 ! anti-symmetric and symmetric terms
   real(wp), external :: Te0_eV, Ti0_eV, density
   real(wp) intpl_FLA_Te0(SG%dim_gq,SG%dim_a)
   real(wp) tmp

   dim_a =SG%dim_a
   dim_gq=SG%dim_gq

   ! call Initialize_MGY_PETSc_Interface_Data(dim_gq,dim_a,ntoroidal)

   ! Tensors on the poloidal plane, and Coeff Tensor including neighboring left and right poloidal plane
   ! Note that there are ghost cells for easier indexing
   allocate(Polar_T(3,3,0:dim_gq,0:dim_a+1))
   allocate(Eps_T(3,3,0:dim_gq,0:dim_a+1,-1:1), CT(3,3,0:dim_gq,0:dim_a+1,-1:1))
   call Determine_GK_PoissonTensors_on_Structured_Grid(0)

   max_nnz=19 ! maximum number of non-zero values per each row
   allocate(tarGid(max_nnz,Lsize), tarVal(max_nnz,Lsize))

   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   tarVal=0.0
   do Lid=1,Lsize
      Gid=mapL2G(Lid)

      qq = mapG2P(Gid)%gqid
      ss = mapG2P(Gid)%aid
      ff = mapG2P(Gid)%gfid

      inv_Jacob=1.0/SG%Jacob_sqf(qq,ss)

      ! 2 (gqid-1) & 3 (gqid+1)
      AntiSym1 = 0.5*(CT(1,2,qq+1,ss,0)-CT(1,2,qq-1,ss,0))
      AntiSym2 = 0.5*(CT(3,2,qq+1,ss,0)-CT(3,2,qq-1,ss,0))
      tarGid(2 ,Lid) = mapP2G(qq-1,ss,ff)
      tarVal(2 ,Lid) = 0.5*(CT(2,2,qq,ss,0)+CT(2,2,qq-1,ss,0)) - min(0.0,AntiSym1) - min(0.0,AntiSym2)
      tarGid(3 ,Lid) = mapP2G(qq+1,ss,ff)
      tarVal(3 ,Lid) = 0.5*(CT(2,2,qq,ss,0)+CT(2,2,qq+1,ss,0)) + max(0.0,AntiSym1) + max(0.0,AntiSym2)

      ! 4 (sid-1) & 5 (sid+1)
      AntiSym1 = 0.5*(CT(2,1,qq,ss+1,0)-CT(2,1,qq,ss-1,0))
      AntiSym2 = 0.5*(CT(3,1,qq,ss+1,0)-CT(3,1,qq,ss-1,0))
      tarGid(4 ,Lid) = mapP2G(qq,ss-1,ff)
      tarVal(4 ,Lid) = 0.5*(CT(1,1,qq,ss,0)+CT(1,1,qq,ss-1,0)) - min(0.0,AntiSym1) - min(0.0,AntiSym2)
      tarGid(5 ,Lid) = mapP2G(qq,ss+1,ff)
      tarVal(5 ,Lid) = 0.5*(CT(1,1,qq,ss,0)+CT(1,1,qq,ss+1,0)) + max(0.0,AntiSym1) + max(0.0,AntiSym2)

      ! 6 (gfid-1) & 7 (gfid+1)
      AntiSym1 = 0.5*(CT(1,3,qq,ss,+1)-CT(1,3,qq,ss,-1))
      AntiSym2 = 0.5*(CT(2,3,qq,ss,+1)-CT(2,3,qq,ss,-1))
      tarGid(6 ,Lid) = mapP2G(qq,ss,ff-1)
      tarVal(6 ,Lid) = 0.5*(CT(3,3,qq,ss,0)+CT(3,3,qq,ss,-1)) - min(0.0,AntiSym1) - min(0.0,AntiSym2)
      tarGid(7 ,Lid) = mapP2G(qq,ss,ff+1)
      tarVal(7 ,Lid) = 0.5*(CT(3,3,qq,ss,0)+CT(3,3,qq,ss,+1)) + max(0.0,AntiSym1) + max(0.0,AntiSym2)




      Sym1 = 0.25*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0))
      ! 8 (sid+1,qid+1) & 9 (sid-1,qid-1)
      tarGid(8 ,Lid) = mapP2G(qq+1,ss+1,ff)
      tarVal(8 ,Lid) = +Sym1
      tarGid(9 ,Lid) = mapP2G(qq-1,ss-1,ff)
      tarVal(9 ,Lid) = +Sym1
      ! 10 (sid+1,qid-1) & 11 (sid-1,qid+1)
      tarGid(10,Lid) = mapP2G(qq-1,ss+1,ff)
      tarVal(10,Lid) = -Sym1
      tarGid(11,Lid) = mapP2G(qq+1,ss-1,ff)
      tarVal(11,Lid) = -Sym1

      Sym1 = 0.25*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0))
      ! 12 (sid+1,gfid+1) & 13 (sid-1,gfid-1)
      tarGid(12,Lid) = mapP2G(qq,ss+1,ff+1)
      tarVal(12,Lid) = +Sym1
      tarGid(13,Lid) = mapP2G(qq,ss-1,ff-1)
      tarVal(13,Lid) = +Sym1
      ! 14 (sid+1,gfid-1) & 15 (sid-1,gfid+1)
      tarGid(14,Lid) = mapP2G(qq,ss+1,ff-1)
      tarVal(14,Lid) = -Sym1
      tarGid(15,Lid) = mapP2G(qq,ss-1,ff+1)
      tarVal(15,Lid) = -Sym1


      Sym1 = 0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))
      ! 16 (gqid+1,gfid+1) & 17  (gqid-1, gfid-1)
      tarGid(16,Lid) = mapP2G(qq+1,ss,ff+1)
      tarVal(16,Lid) = +Sym1
      tarGid(17,Lid) = mapP2G(qq-1,ss,ff-1)
      tarVal(17,Lid) = +Sym1
      ! 18 (gqid+1,gfid-1) & 18 (gqid-1, gfid+1)
      tarGid(18,Lid) = mapP2G(qq+1,ss,ff-1)
      tarVal(18,Lid) = -Sym1
      tarGid(19,Lid) = mapP2G(qq-1,ss,ff+1)
      tarVal(19,Lid) = -Sym1


      ! 1 (gqid,sid,gfid)
      tarGid(1,Lid)=mapP2G(qq,ss,ff)
      tarVal(1,Lid)=-sum(tarVal(2:19,Lid))

      ! IMPORTANT!! Multiply -inv_Jacob to the above coefficients,
      ! so make sure that target coefficient values represent (-nabla*(eps_T*nabla*phi))
      tarVal(:,Lid) = -1.0*inv_Jacob*tarVal(:,Lid)

   enddo
   call Apply_BC_to_GK_Poisson_LHS_Coeffs ! Apply Boundary Condition to Amt Coeffs
   call Add_Adiabatic_Term_to_GK_Poisson_LHS



   ! To calculate numbers of non-zero elements for diagonal(d_nnz) & off-diagonal(o_nnz) blocks
   allocate(d_nnz(Lsize), o_nnz(Lsize))
   d_nnz=0; o_nnz=0;
   do Lid=1,Lsize
      do tnum=1,max_nnz
         if(tarGid(tnum,Lid)>=0) then
            if(tarGid(tnum,Lid)<Gid_start .or. tarGid(tnum,Lid)>Gid_end) then
               o_nnz(Lid)=o_nnz(Lid)+1
            else
               d_nnz(Lid)=d_nnz(Lid)+1
            endif
         endif
      enddo
   enddo


   ! Construct Matrix (PMat)
   call MatCreate(PETSC_COMM_WORLD,PMat,ierr); CHKERRQ(ierr)
   call MatSetSizes(PMat,Lsize,Lsize,Gsize,Gsize,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(PMat,ierr); CHKERRQ(ierr)
   ! Preallocate the memory using o_nnz and d_nnz
   call MatMPIAIJSetPreallocation(PMat,max_nnz,d_nnz,max_nnz,o_nnz,ierr); CHKERRQ(ierr)
   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(PMat,1,Gid,max_nnz,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "PMat Assembly takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call MatView(PMat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)


   call VecCreate(PETSC_COMM_WORLD,bVec,ierr); CHKERRQ(ierr)
   call VecSetSizes(bVec,Lsize,Gsize,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(bVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,xVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,ResVec,ierr); CHKERRQ(ierr)

   if(Flag_Def_dne_adiab==3) then

      ! Set PMat value as the first term
      call MatDuplicate(PMat,MAT_COPY_VALUES,Jmat,ierr); CHKERRQ(ierr)

      ! Add another term from the nonlinear adiabatic response
      do Lid=1,Lsize
         Gid=mapL2G(Lid)
         qq = mapG2P(Gid)%gqid
         ss = mapG2P(Gid)%aid
         tmp = 1e6*SG%n0_1D(ss)*(1.602176634e-19)/SG%Te0_eV_1D(ss)*exp(SG%phi(qq,ss,1)/SG%Te0_eV_1D(ss))
         call MatSetValues(Jmat,1,Gid,1,Gid,tmp,ADD_VALUES,ierr); CHKERRQ(ierr)
      enddo

      ! Assembly Jacobian Matrix
      tt0=MPI_WTIME()
      call MatAssemblyBegin(Jmat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
      call MatAssemblyEnd(Jmat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "JMat Assembly takes ",tt1-tt0," secs"
      call flush(stdout)

      !pre_Jmat
      call MatDuplicate(Jmat,MAT_COPY_VALUES,pre_Jmat,ierr); CHKERRQ(ierr)

      ! call MatView(Jmat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

      call VecDuplicate(bVec,srcVec,ierr); CHKERRQ(ierr)
      call VecDuplicate(bVec,yVec,ierr); CHKERRQ(ierr)
   endif


   if(Flag_Def_dne_adiab==3) then
      ! Construct KSP
      tt0=MPI_WTIME()
      call KSPCreate(PETSC_COMM_WORLD,ksp_P,ierr); CHKERRQ(ierr)
      ! call KSPSetOperators(ksp_P,Jmat,Jmat,ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp_P,Jmat,pre_Jmat,ierr); CHKERRQ(ierr)
      ! call KSPSetTolerances(ksp_P,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
      call KSPSetFromOptions(ksp_P,ierr); CHKERRQ(ierr)
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
      call flush(stdout)
   else
      ! Construct KSP
      tt0=MPI_WTIME()
      call KSPCreate(PETSC_COMM_WORLD,ksp_P,ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp_P,PMat,PMat,ierr); CHKERRQ(ierr)
      ! call KSPSetTolerances(ksp_P,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
      call KSPSetFromOptions(ksp_P,ierr); CHKERRQ(ierr)
      tt1=MPI_WTIME()
      if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
      call flush(stdout)

   endif

end subroutine Construct_GK_Poisson3D_Solver_with_Extremal_Method


subroutine Update_GK_Poisson_LHS_Coeffs(In_Flag_WhichLHS)

   use Structured_Grid
   use particle_decomp
   use global_parameters
   use field_array, only:FLA_Te0
   implicit none

   integer,intent(in),optional :: In_Flag_WhichLHS
   integer :: Flag_WhichLHS
   PetscInt Lid,Gid,qq,ss,ff,tnum
   real tt0, tt1
   PetscReal inv_Jacob
   PetscErrorCode ierr
   integer :: dim_gq, dim_a, i
   PetscReal AntiSym1, AntiSym2, Sym1 ! anti-symmetric and symmetric terms
   real(wp), external :: Te0_eV, Ti0_eV, density
   real(wp) intpl_FLA_Te0(SG%dim_gq,SG%dim_a)


   Flag_WhichLHS=0 !Default value=0 (Full Eps_T); =1(Metric_T only); =2(Polar_T only)
   if(present(In_Flag_WhichLHS)) Flag_WhichLHS = In_Flag_WhichLHS

   if(Flag_WhichLHS==0 .or. Flag_WhichLHS==1 .or. Flag_WhichLHS==2) then
   else
      write(stdout,*) "Error!! In Update_GK_Poisson_LHS_Coeffs, Flag_WhichLHS=",Flag_WhichLHS
      call flush(stdout)
      stop
   endif

   dim_a =SG%dim_a
   dim_gq=SG%dim_gq

   ! Calculate dn
   if(Flag_Weight_dfh==0) then
      SG%dne_SM = SG%dne(:,:,1)
      SG%dni_SM = SG%dni(:,:,1)
   elseif(Flag_Weight_dfh==1) then
      do i=1,SG%dim_s
         SG%dne_SM(:,i) = SG%ne_by_fh(:,i,1) -SG%n0_1D(i)
         SG%dni_SM(:,i) = SG%ni_by_ff(:,i,1) -SG%n0_1D(i)
      enddo
   elseif(Flag_Weight_dfh==2 .or. Flag_Weight_dfh==3) then
      ! dn = Hybrid_Density - n0
      ! where Hybrid_Density = hyb_ratio*(full_f or full_h) + (1-hyb_ratio)*(delta_f or delta_h)
      do i=1,SG%dim_s
         SG%dne_SM(:,i) = SG%hyb_ratio_e(:,i)*SG%ne_by_fh(:,i,1)+(1.0-SG%hyb_ratio_e(:,i))*(SG%dne(:,i,1)+SG%n0_1D(i)) -SG%n0_1D(i)
         SG%dni_SM(:,i) = SG%hyb_ratio_i(:,i)*SG%ni_by_ff(:,i,1)+(1.0-SG%hyb_ratio_i(:,i))*(SG%dni(:,i,1)+SG%n0_1D(i)) -SG%n0_1D(i)
      enddo
   endif

   ! Then, smooth them
   call Smooth_SG_Data_2D([2,0],SG%dne_SM,SG%Cell_dV)
   call Smooth_SG_Data_2D([2,0],SG%dni_SM,SG%Cell_dV)


   call Determine_GK_PoissonTensors_on_Structured_Grid(Flag_WhichLHS)

   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   tarVal=0.0
   do Lid=1,Lsize
      Gid=mapL2G(Lid)

      qq = mapG2P(Gid)%gqid
      ss = mapG2P(Gid)%aid
      ff = mapG2P(Gid)%gfid

      inv_Jacob=1.0/SG%Jacob_sqf(qq,ss)

      ! 2 (gqid-1) & 3 (gqid+1)
      AntiSym1 = 0.5*(CT(1,2,qq+1,ss,0)-CT(1,2,qq-1,ss,0))
      AntiSym2 = 0.5*(CT(3,2,qq+1,ss,0)-CT(3,2,qq-1,ss,0))
      tarGid(2 ,Lid) = mapP2G(qq-1,ss,ff)
      tarVal(2 ,Lid) = 0.5*(CT(2,2,qq,ss,0)+CT(2,2,qq-1,ss,0)) - min(0.0,AntiSym1) - min(0.0,AntiSym2)
      tarGid(3 ,Lid) = mapP2G(qq+1,ss,ff)
      tarVal(3 ,Lid) = 0.5*(CT(2,2,qq,ss,0)+CT(2,2,qq+1,ss,0)) + max(0.0,AntiSym1) + max(0.0,AntiSym2)

      ! 4 (sid-1) & 5 (sid+1)
      AntiSym1 = 0.5*(CT(2,1,qq,ss+1,0)-CT(2,1,qq,ss-1,0))
      AntiSym2 = 0.5*(CT(3,1,qq,ss+1,0)-CT(3,1,qq,ss-1,0))
      tarGid(4 ,Lid) = mapP2G(qq,ss-1,ff)
      tarVal(4 ,Lid) = 0.5*(CT(1,1,qq,ss,0)+CT(1,1,qq,ss-1,0)) - min(0.0,AntiSym1) - min(0.0,AntiSym2)
      tarGid(5 ,Lid) = mapP2G(qq,ss+1,ff)
      tarVal(5 ,Lid) = 0.5*(CT(1,1,qq,ss,0)+CT(1,1,qq,ss+1,0)) + max(0.0,AntiSym1) + max(0.0,AntiSym2)

      ! 6 (gfid-1) & 7 (gfid+1)
      AntiSym1 = 0.5*(CT(1,3,qq,ss,+1)-CT(1,3,qq,ss,-1))
      AntiSym2 = 0.5*(CT(2,3,qq,ss,+1)-CT(2,3,qq,ss,-1))
      tarGid(6 ,Lid) = mapP2G(qq,ss,ff-1)
      tarVal(6 ,Lid) = 0.5*(CT(3,3,qq,ss,0)+CT(3,3,qq,ss,-1)) - min(0.0,AntiSym1) - min(0.0,AntiSym2)
      tarGid(7 ,Lid) = mapP2G(qq,ss,ff+1)
      tarVal(7 ,Lid) = 0.5*(CT(3,3,qq,ss,0)+CT(3,3,qq,ss,+1)) + max(0.0,AntiSym1) + max(0.0,AntiSym2)


      Sym1 = 0.25*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0))
      ! 8 (sid+1,qid+1) & 9 (sid-1,qid-1)
      tarGid(8 ,Lid) = mapP2G(qq+1,ss+1,ff)
      tarVal(8 ,Lid) = +Sym1
      tarGid(9 ,Lid) = mapP2G(qq-1,ss-1,ff)
      tarVal(9 ,Lid) = +Sym1
      ! 10 (sid+1,qid-1) & 11 (sid-1,qid+1)
      tarGid(10,Lid) = mapP2G(qq-1,ss+1,ff)
      tarVal(10,Lid) = -Sym1
      tarGid(11,Lid) = mapP2G(qq+1,ss-1,ff)
      tarVal(11,Lid) = -Sym1

      Sym1 = 0.25*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0))
      ! 12 (sid+1,gfid+1) & 13 (sid-1,gfid-1)
      tarGid(12,Lid) = mapP2G(qq,ss+1,ff+1)
      tarVal(12,Lid) = +Sym1
      tarGid(13,Lid) = mapP2G(qq,ss-1,ff-1)
      tarVal(13,Lid) = +Sym1
      ! 14 (sid+1,gfid-1) & 15 (sid-1,gfid+1)
      tarGid(14,Lid) = mapP2G(qq,ss+1,ff-1)
      tarVal(14,Lid) = -Sym1
      tarGid(15,Lid) = mapP2G(qq,ss-1,ff+1)
      tarVal(15,Lid) = -Sym1


      Sym1 = 0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))
      ! 16 (gqid+1,gfid+1) & 17  (gqid-1, gfid-1)
      tarGid(16,Lid) = mapP2G(qq+1,ss,ff+1)
      tarVal(16,Lid) = +Sym1
      tarGid(17,Lid) = mapP2G(qq-1,ss,ff-1)
      tarVal(17,Lid) = +Sym1
      ! 18 (gqid+1,gfid-1) & 18 (gqid-1, gfid+1)
      tarGid(18,Lid) = mapP2G(qq+1,ss,ff-1)
      tarVal(18,Lid) = -Sym1
      tarGid(19,Lid) = mapP2G(qq-1,ss,ff+1)
      tarVal(19,Lid) = -Sym1


      ! 1 (gqid,sid,gfid)
      tarGid(1,Lid)=mapP2G(qq,ss,ff)
      tarVal(1,Lid)=-sum(tarVal(2:19,Lid))

      ! IMPORTANT!! Multiply -inv_Jacob to the above coefficients,
      ! so make sure that target coefficient values represent (-nabla*(eps_T*nabla*phi))
      tarVal(:,Lid) = -1.0*inv_Jacob*tarVal(:,Lid)

   enddo
   call Apply_BC_to_GK_Poisson_LHS_Coeffs ! Apply Boundary Condition to Amt Coeffs


   if(Flag_WhichLHS==0) then
      call Add_Adiabatic_Term_to_GK_Poisson_LHS
   endif


   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(PMat,1,Gid,max_nnz,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(PMat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "PMat Assembly takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call MatView(PMat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)


end subroutine Update_GK_Poisson_LHS_Coeffs

subroutine Determine_GK_PoissonTensors_on_Structured_Grid(In_Flag_WhichTerm)
   use precision, only:wp
   use global_parameters, only:stdout,mype,ifelectron,mzetamax,msm0,amasse,amassi
   use Structured_Grid
   use particle_decomp
   implicit none

   integer, intent(in), optional :: In_Flag_WhichTerm
   integer Flag_WhichTerm
   integer sid,gqid,i,j,ierr
   real(wp) invBsq,tmp ! 1/(B^2)
   real(wp), external :: density
   real(wp) inv_du(3)
   integer istatus(MPI_STATUS_SIZE), icount, ierror
   real(wp) ne_tot(SG%dim_gq,SG%dim_s), ni_tot(SG%dim_gq,SG%dim_s)

   Flag_WhichTerm=2 !Default value=0 (Full Eps_T); =1(Metric_T only); =2(Polar_T only)
   if(present(In_Flag_WhichTerm)) Flag_WhichTerm = In_Flag_WhichTerm

   if(Flag_Def_dne_adiab==1) then
      do sid=1,SG%dim_s
         do gqid=1,SG%dim_gq
            ! density unit is (#/m^3)
            ! ni_tot(gqid,sid) = 1e6*(SG%n0_1D(sid) + SG%dni(gqid,sid,1))
            ! ne_tot(gqid,sid) = 1e6*(SG%n0_1D(sid) + SG%dne(gqid,sid,1) + int(ifelectron/2)*SG%phi(gqid,sid,1)/(SG%Te0_eV_1D(sid))*SG%n0_1D(sid))
            ni_tot(gqid,sid) = 1e6*(SG%n0_1D(sid) + SG%dni_SM(gqid,sid))
            ne_tot(gqid,sid) = 1e6*(SG%n0_1D(sid) + SG%dne_SM(gqid,sid) + int(ifelectron/2)*SG%phi(gqid,sid,1)/(SG%Te0_eV_1D(sid))*SG%n0_1D(sid))
         enddo
      enddo
   elseif(Flag_Def_dne_adiab==3) then
      ! ne_tot=n0*exp(e*Phi/Te0) + dne_NA
      do sid=1,SG%dim_s
         do gqid=1,SG%dim_gq
            ! density unit is (#/m^3)
            ni_tot(gqid,sid) = 1e6*(SG%n0_1D(sid) + SG%dni_SM(gqid,sid))
            ne_tot(gqid,sid) = 1e6*(SG%dne_SM(gqid,sid) + SG%n0_1D(sid)*exp(int(ifelectron/2)*SG%phi(gqid,sid,1)/(SG%Te0_eV_1D(sid))))
         enddo
      enddo

   else
      write(stdout,*) "Error! Flag_Def_dne_adiab=", Flag_Def_dne_adiab, " in Determine_GK_PoissonTensor"
      call flush(stdout)
      stop
   endif

   ! Density must be positive so ignore negative values because of numerical noises
   ne_tot = max(0.0,ne_tot)
   ni_tot = max(0.0,ni_tot)

   do sid=1,SG%dim_a
      ! Assuming the static equilibrium density for ions and electrons
      do gqid=1,SG%dim_gq
         invBsq=1.0/(SG%BB(gqid,sid)**2.0) ! unit (T^2)
         ! write(stdout,"(3(A,ES16.8))") "n0(sid)=",n0,", B=",SG%BB(gqid,sid), ", Polar(sid)=", 1e-3*(n0*amasse + n0*amassi)*invBsq
         do j=1,3
            do i=1,3
               ! Note tsince mgy_petsc coefficients are using SI units
               ! Following (1e-3) is a converting coefficient for mass (g->kg)
               tmp=1e-3*(ne_tot(gqid,sid)*amasse + ni_tot(gqid,sid)*amassi)*invBsq
               Polar_T(i,j,gqid,sid)= tmp*(Metric_T(i,j,gqid,sid)- invBsq*(SG%B0_cntv3_sqf(i,gqid,sid)*SG%B0_cntv3_sqf(j,gqid,sid)))
            enddo
         enddo
      enddo
   enddo
   Polar_T(:,:,SG%dim_gq,:)=Polar_T(:,:,1,:)
   Polar_T(:,:,0,:)        =Polar_T(:,:,SG%dim_gq-1,:)
   Polar_T(:,:,:,0)        =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   Polar_T(:,:,:,SG%dim_a+1)  =-100 ! Assign Negative initial values which will be ignored in MatSetValues

   if(Flag_WhichTerm==0) then
      ! Construct Epsilon Tensor ( eps=eps_0*Metric_T + Sigma_over_s(ns*ms/B^2 * (Metric_T - B^i*B^j/B^2)) )
      Eps_T(:,:,:,:,0)=8.8541878128e-12*Metric_T+Polar_T
   elseif(Flag_WhichTerm==1) then
      ! Eps=eps_0*Metric_T
      Eps_T(:,:,:,:,0)=8.8541878128e-12*Metric_T
   elseif(Flag_WhichTerm==2) then
      ! Eps=Polar_T
      Eps_T(:,:,:,:,0)=Polar_T
   else
      write(stdout,*) "Error!! Flag_WhichTerm=", Flag_WhichTerm
      call flush(stdout)
      stop
   endif

   ! Communicate Eps_T with neighboring left and right planes.
   ! So, Eps_T(:,:,:,:,tid) ==>  tid=-1(left plane), tid=0(my plane), tid=1(right plane)
   ! Eps_T(3,3,0:dim_gq,0:dim_a+1,-1:1)
   icount=3*3*(SG%dim_gq+1)*(SG%dim_a+2)
   call MPI_SENDRECV(Eps_T(:,:,:,:,0),icount,mpi_Rsize,left_pe,0,Eps_T(:,:,:,:,1),icount,mpi_Rsize,right_pe,0,toroidal_comm,istatus,ierror)
   call MPI_SENDRECV(Eps_T(:,:,:,:,0),icount,mpi_Rsize,right_pe,0,Eps_T(:,:,:,:,-1),icount,mpi_Rsize,left_pe,0,toroidal_comm,istatus,ierror)



   inv_du(1) = SG%inv_ds
   inv_du(2) = SG%inv_dgq
   inv_du(3) = 1.0/(SG%zetamax-SG%zetamin)

   ! Determine Coefficient Tensor (CT = J*Eps_T(i,j)/(du^i*du^j))
   do sid=1,SG%dim_a
      do gqid=1,SG%dim_gq
         do j=1,3
            do i=1,3
               CT(i,j,gqid,sid,:) = SG%Jacob_sqf(gqid,sid)*Eps_T(i,j,gqid,sid,:)*inv_du(i)*inv_du(j)
            enddo
         enddo
      enddo
   enddo
   CT(:,:,SG%dim_gq,:,:)  =CT(:,:,1,:,:)
   CT(:,:,0,:,:)          =CT(:,:,SG%dim_gq-1,:,:)
   CT(:,:,:,0,:)          =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   CT(:,:,:,SG%dim_a+1,:) =-100 ! Assign Negative initial values which will be ignored in MatSetValues


end subroutine Determine_GK_PoissonTensors_on_Structured_Grid


subroutine Add_Adiabatic_Term_to_GK_Poisson_LHS

   use precision
   use Structured_Grid
   use field_array

   implicit none

   integer Lid,Gid,qq,ss
   real(wp), external :: density, Te0_eV
   real(wp) intpl_FLA_Te0(SG%dim_gq,SG%dim_a)
   real(wp) AD_coeff


   !In the case of delta-h electron simulation,
   !add the adiabatic electron response term in the LHS (PMat)
   ! NOTE: The following line must locate after "call Apply_BC_to_GK_Poisson_LHS_Coeffs"
   ! Note that, in the case of Flag_Def_dne_adiab==3, there is no add to tarVal, so tarVal simply represents -nabla*eps*nabla differential operator
   if(ifelectron==2) then
      if(Flag_Def_dne_adiab==0 .or. Flag_Def_dne_adiab==1) then
         do Lid=1,Lsize
            Gid= mapL2G(Lid)
            qq = mapG2P(Gid)%gqid
            ss = mapG2P(Gid)%aid

            AD_coeff= (1.602176634e-19)*1e6*density(SG%a_1D(ss),SG%gq_1D(qq))/Te0_eV(SG%a_1D(ss),SG%gq_1D(qq))


            tarVal(1,Lid)=tarVal(1,Lid) + AD_coeff! 1 (Center)

            ! Adiabatic Density is not only determined by the phi on that node but also affected by neighboring phi
            ! (See a note "The Poisson equation regarding the electron adiabatic density")
            ! Determine tarVal corresponding to the integral coefficients of the shape function
            if(ss==1) then
               ! ! Left Half case
               ! AD_coeff=AD_coeff/18.0
               ! tarVal(1,Lid)=tarVal(1,Lid) + 8*AD_coeff ! 1 (Center)
               ! tarVal(2,Lid)=tarVal(2,Lid) + 2*AD_coeff ! 2 (qq-1)
               ! tarVal(3,Lid)=tarVal(3,Lid) + 2*AD_coeff ! 3 (qq+1)
               ! tarVal(5,Lid)=tarVal(5,Lid) + 4*AD_coeff ! 5 (ss+1)
               ! tarVal(8,Lid)=tarVal(8,Lid) + 1*AD_coeff ! 8 (ss+1,qq+1)
               ! tarVal(10,Lid)=tarVal(10,Lid) + 1*AD_coeff ! 10 (ss+1,qq-1)
            elseif(ss==SG%dim_s) then
               ! ! Right Half case
               ! AD_coeff=AD_coeff/18.0
               ! tarVal(1,Lid)=tarVal(1,Lid) + 8*AD_coeff ! 1 (Center)
               ! tarVal(2,Lid)=tarVal(2,Lid) + 2*AD_coeff ! 2 (qq-1)
               ! tarVal(3,Lid)=tarVal(3,Lid) + 2*AD_coeff ! 3 (qq+1)
               ! tarVal(4,Lid)=tarVal(4,Lid) + 4*AD_coeff ! 4 (ss-1)
               ! tarVal(9,Lid)=tarVal(9,Lid) + 1*AD_coeff ! 9 (ss-1,qq-1)
               ! tarVal(11,Lid)=tarVal(11,Lid) + 1*AD_coeff ! 11 (ss-1,qq+1)
            else
               ! ! Full case
               ! AD_coeff=AD_coeff/36.0
               ! tarVal(1,Lid)=tarVal(1,Lid) + 16*AD_coeff! 1 (Center)
               ! tarVal(2,Lid)=tarVal(2,Lid) + 4*AD_coeff ! 2 (qq-1)
               ! tarVal(3,Lid)=tarVal(3,Lid) + 4*AD_coeff ! 3 (qq+1)
               ! tarVal(4,Lid)=tarVal(4,Lid) + 4*AD_coeff ! 4 (ss-1)
               ! tarVal(5,Lid)=tarVal(5,Lid) + 4*AD_coeff ! 5 (ss+1)
               ! tarVal(8,Lid)=tarVal(8,Lid) + 1*AD_coeff ! 8 (ss+1,qq+1)
               ! tarVal(9,Lid)=tarVal(9,Lid) + 1*AD_coeff ! 9 (ss-1,qq-1)
               ! tarVal(10,Lid)=tarVal(10,Lid) + 1*AD_coeff ! 10 (ss+1,qq-1)
               ! tarVal(11,Lid)=tarVal(11,Lid) + 1*AD_coeff ! 11 (ss-1,qq+1)
            endif
         enddo
      elseif(Flag_Def_dne_adiab==2) then
         ! use FLA_Te0 to the definition of the electron adiabatic density (dne_adiab = n0*e*full_phi/FLA_Te0)
         call Interpolate_OriGridData_into_SG(FLA_Te0(1,:), intpl_FLA_Te0)
         do Lid=1,Lsize
            Gid= mapL2G(Lid)
            qq = mapG2P(Gid)%gqid
            ss = mapG2P(Gid)%aid
            tarVal(1,Lid)=tarVal(1,Lid) + (1.602176634e-19)*1e6*density(SG%a_1D(ss),SG%gq_1D(qq))/(1e3*Tem_00*intpl_FLA_Te0(qq,ss))
         enddo
      endif

   elseif(ifelectron==1 .and. Flag_Poisson3D_on_SG==2) then
      ! To mimic GTS Way
      if(adiab > 0.0_wp) then
         do Lid=1,Lsize
            Gid= mapL2G(Lid)
            qq = mapG2P(Gid)%gqid
            ss = mapG2P(Gid)%aid
            tarVal(1,Lid)=tarVal(1,Lid) + adiab*(1.602176634e-19)*1e6*density(SG%a_1D(ss),SG%gq_1D(qq))/Te0_eV(SG%a_1D(ss),SG%gq_1D(qq))
         enddo
      else
         do Lid=1,Lsize
            Gid= mapL2G(Lid)
            qq = mapG2P(Gid)%gqid
            ss = mapG2P(Gid)%aid
            tarVal(1,Lid)=tarVal(1,Lid) + 0.01/(Cs*Cs)/(1e3*Tem_00)*(1.602176634e-19)*1e6*density(SG%a_1D(ss),SG%gq_1D(qq))
         enddo
      endif
   endif


end subroutine Add_Adiabatic_Term_to_GK_Poisson_LHS


subroutine Apply_BC_to_GK_Poisson_LHS_Coeffs
   use Structured_Grid
   use global_parameters, only:mzetamax,stdout,mype
   implicit none

   PetscInt Lid,Gid,qq,ss,ff,j
   integer :: dim_gq, dim_a
   PetscReal inv_Jacob, Kj_coeff(3), inv_du(3)

   dim_a =SG%dim_a
   dim_gq=SG%dim_gq

   inv_du(1) = SG%inv_ds
   inv_du(2) = SG%inv_dgq
   inv_du(3) = 1.0/(SG%zetamax-SG%zetamin)

   ! Apply Neumann (+Dirichlet) Condition at inner and outer boundary
   if(BC0_type==0) then !Dirichlet condition
      ss=1
      do ff=1,mzetamax
         do qq=1,dim_gq-1

            !Dirichlet condition
            Lid=mapP2L(qq,ss,ff)
            if(Lid>0) then
               tarGid(:,Lid) = -100
               tarVal(:,Lid) = 0.0

               tarGid(1,Lid) = mapP2G(qq,ss,ff)
               tarVal(1,Lid) = 1.0 !phi @ (ss=1)
            endif
         enddo
      enddo
   elseif(BC0_type==1) then
      ss=1
      do ff=1,mzetamax
         do qq=1,dim_gq-1
            ! Neumann condition (abs(E)=0 at inner)
            Lid=mapP2L(qq,ss,ff)
            if(Lid>0) then

               inv_Jacob=1.0/SG%Jacob_sqf(qq,ss)

               do j=1,3
                  Kj_coeff(j) = 0.5*inv_du(1)*(-3*SG%Jacob_sqf(qq,ss)*Eps_T(1,j,qq,ss,0) +4*SG%Jacob_sqf(qq,ss+1)*Eps_T(1,j,qq,ss+1,0) -SG%Jacob_sqf(qq,ss+2)*Eps_T(1,j,qq,ss+2,0)) &
                              +0.5*inv_du(2)*(SG%Jacob_sqf(qq+1,ss)*Eps_T(2,j,qq+1,ss,0) -SG%Jacob_sqf(qq-1,ss)*Eps_T(2,j,qq-1,ss,0)) &
                              +0.5*inv_du(3)*(SG%Jacob_sqf(qq,ss)*Eps_T(3,j,qq,ss,+1) -SG%Jacob_sqf(qq,ss)*Eps_T(3,j,qq,ss,-1));
               enddo

               tarGid(:, Lid) = -100
               tarVal(:, Lid) = 0.0

               tarGid(2, Lid) = mapP2G(qq-1,ss,ff)
               tarVal(2, Lid) = CT(2,2,qq,ss,0) + 0.5*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0)) - 0.5*Kj_coeff(2)*inv_du(2)
               tarGid(3, Lid) = mapP2G(qq+1,ss,ff)
               tarVal(3, Lid) = CT(2,2,qq,ss,0) - 0.5*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0)) + 0.5*Kj_coeff(2)*inv_du(2)

               tarGid(4, Lid) = mapP2G(qq,ss+1,ff)
               tarVal(4, Lid) = 2.0*CT(1,1,qq,ss,0)

               tarGid(5, Lid) = mapP2G(qq,ss,ff-1)
               tarVal(5, Lid) = CT(3,3,qq,ss,0) + 0.5*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0)) - 0.5*Kj_coeff(3)*inv_du(3)
               tarGid(6, Lid) = mapP2G(qq,ss,ff+1)
               tarVal(6, Lid) = CT(3,3,qq,ss,0) - 0.5*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0)) + 0.5*Kj_coeff(3)*inv_du(3)

               tarGid(7, Lid) = mapP2G(qq+1,ss+1,ff)
               tarVal(7, Lid) = +0.5*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0))
               tarGid(8, Lid) = mapP2G(qq-1,ss+1,ff)
               tarVal(8, Lid) = -0.5*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0))

               tarGid(9, Lid) = mapP2G(qq,ss+1,ff+1)
               tarVal(9, Lid) = +0.5*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0))
               tarGid(10,Lid) = mapP2G(qq,ss+1,ff-1)
               tarVal(10,Lid) = -0.5*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0))

               tarGid(11,Lid) = mapP2G(qq+1,ss,ff+1)
               tarVal(11,Lid) = +0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))
               tarGid(12,Lid) = mapP2G(qq-1,ss,ff-1)
               tarVal(12,Lid) = +0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))

               tarGid(13,Lid) = mapP2G(qq+1,ss,ff-1)
               tarVal(13,Lid) = -0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))
               tarGid(14,Lid) = mapP2G(qq-1,ss,ff+1)
               tarVal(14,Lid) = -0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))

               tarGid(1, Lid) = mapP2G(qq,ss,ff)
               tarVal(1, Lid) = -sum(tarVal(2:14,Lid))


               ! IMPORTANT!! Multiply -inv_Jacob to the above coefficients,
               ! so make sure that target coefficient values represent (-nabla*(eps_T*nabla*phi))
               tarVal(:,Lid) = -1.0*inv_Jacob*tarVal(:,Lid)

            endif
         enddo
      enddo
   endif


   if(BC1_type==0) then
      ss=dim_a
      do ff=1,mzetamax
         do qq=1,dim_gq-1
            Lid=mapP2L(qq,ss,ff)
            if(Lid>0) then
               tarGid(:, Lid) = -100
               tarVal(:, Lid) = 0.0

               tarGid(1, Lid) = mapP2G(qq,ss,ff)
               tarVal(1, Lid) = 1.0 !phi(ss=dim_a)
            endif
         enddo
      enddo
   elseif(BC1_type==1) then
      ss=dim_a
      inv_du(1)=-1.0*SG%inv_ds ! since delta_s is now negative value for outer boundary case
      do ff=1,mzetamax
         do qq=1,dim_gq-1
            ! Neumann condition (abs(E)=0 at inner)
            Lid=mapP2L(qq,ss,ff)
            if(Lid>0) then

               inv_Jacob=1.0/SG%Jacob_sqf(qq,ss)

               do j=1,3
                  Kj_coeff(j) = 0.5*inv_du(1)*(-3*SG%Jacob_sqf(qq,ss)*Eps_T(1,j,qq,ss,0) +4*SG%Jacob_sqf(qq,ss-1)*Eps_T(1,j,qq,ss-1,0) -SG%Jacob_sqf(qq,ss-2)*Eps_T(1,j,qq,ss-2,0)) &
                              +0.5*inv_du(2)*(SG%Jacob_sqf(qq+1,ss)*Eps_T(2,j,qq+1,ss,0) -SG%Jacob_sqf(qq-1,ss)*Eps_T(2,j,qq-1,ss,0)) &
                              +0.5*inv_du(3)*(SG%Jacob_sqf(qq,ss)*Eps_T(3,j,qq,ss,+1) -SG%Jacob_sqf(qq,ss)*Eps_T(3,j,qq,ss,-1));
               enddo

               tarGid(:, Lid) = -100
               tarVal(:, Lid) = 0.0

               tarGid(2, Lid) = mapP2G(qq-1,ss,ff)
               tarVal(2, Lid) = CT(2,2,qq,ss,0) + 0.5*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0)) - 0.5*Kj_coeff(2)*inv_du(2)
               tarGid(3, Lid) = mapP2G(qq+1,ss,ff)
               tarVal(3, Lid) = CT(2,2,qq,ss,0) - 0.5*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0)) + 0.5*Kj_coeff(2)*inv_du(2)

               tarGid(4, Lid) = mapP2G(qq,ss-1,ff)
               tarVal(4, Lid) = 2.0*CT(1,1,qq,ss,0)

               tarGid(5, Lid) = mapP2G(qq,ss,ff-1)
               tarVal(5, Lid) = CT(3,3,qq,ss,0) + 0.5*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0)) - 0.5*Kj_coeff(3)*inv_du(3)
               tarGid(6, Lid) = mapP2G(qq,ss,ff+1)
               tarVal(6, Lid) = CT(3,3,qq,ss,0) - 0.5*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0)) + 0.5*Kj_coeff(3)*inv_du(3)

               tarGid(7, Lid) = mapP2G(qq+1,ss-1,ff)
               tarVal(7, Lid) = +0.5*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0))
               tarGid(8, Lid) = mapP2G(qq-1,ss-1,ff)
               tarVal(8, Lid) = -0.5*(CT(1,2,qq,ss,0)+CT(2,1,qq,ss,0))

               tarGid(9, Lid) = mapP2G(qq,ss-1,ff+1)
               tarVal(9, Lid) = +0.5*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0))
               tarGid(10,Lid) = mapP2G(qq,ss-1,ff-1)
               tarVal(10,Lid) = -0.5*(CT(1,3,qq,ss,0)+CT(3,1,qq,ss,0))

               tarGid(11,Lid) = mapP2G(qq+1,ss,ff+1)
               tarVal(11,Lid) = +0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))
               tarGid(12,Lid) = mapP2G(qq-1,ss,ff-1)
               tarVal(12,Lid) = +0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))

               tarGid(13,Lid) = mapP2G(qq+1,ss,ff-1)
               tarVal(13,Lid) = -0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))
               tarGid(14,Lid) = mapP2G(qq-1,ss,ff+1)
               tarVal(14,Lid) = -0.25*(CT(2,3,qq,ss,0)+CT(3,2,qq,ss,0))

               tarGid(1, Lid) = mapP2G(qq,ss,ff)
               tarVal(1, Lid) = -sum(tarVal(2:14,Lid))


               ! IMPORTANT!! Multiply -inv_Jacob to the above coefficients,
               ! so make sure that target coefficient values represent (-nabla*(eps_T*nabla*phi))
               tarVal(:,Lid) = -1.0*inv_Jacob*tarVal(:,Lid)

            endif
         enddo
      enddo
   endif


   if(BC0_type==1 .and. BC1_type==1) then
      ! This is purely Neumann conditions.. so the solution is not unique
      ! We need a reference point (Dirichlet condition)
      Lid=mapP2L(1,1,1) ! This is a reference point
      if(Lid>0) then
         tarGid(:, Lid) = -100
         tarVal(:, Lid) = 0.0

         tarGid(1, Lid) = mapL2G(Lid)
         tarVal(1, Lid) = 1.0
      endif
   endif

end subroutine Apply_BC_to_GK_Poisson_LHS_Coeffs

subroutine Apply_BC_to_GK_Poisson_RHS_Coeffs
   use Structured_Grid
   use global_parameters, only:mzetamax,stdout,mype
   implicit none

   PetscInt Lid,Gid,qq,ss,ff,j
   integer :: dim_gq, dim_a
   PetscReal Kss_coeff, inv_du(3), val
   PetscErrorCode ierr


   dim_a =SG%dim_a
   dim_gq=SG%dim_gq

   inv_du(1) = SG%inv_ds
   inv_du(2) = SG%inv_dgq
   inv_du(3) = 1.0/(SG%zetamax-SG%zetamin)


   if(BC0_type==0) then ! Dirichlet condition (phi at inner boundary)
      ss=1
      do ff=1,mzetamax
         do qq=1,dim_gq-1
            Lid=mapP2L(qq,ss,ff)
            if(Lid>0) then
               Gid=mapL2G(Lid)
               val=BC0_value
               call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
            endif
         enddo
      enddo
   elseif(BC0_type==1) then ! Neumann condition (Er value at inner boundary)
      ss=1
      do ff=1,mzetamax
         do qq=1,dim_gq-1
            Lid=mapP2L(qq,ss,ff)
            if(Lid>0) then
               Gid=mapL2G(Lid)
               Kss_coeff = 0.5*inv_du(1)*(-3*SG%Jacob_sqf(qq,ss)*Eps_T(1,1,qq,ss,0) +4*SG%Jacob_sqf(qq,ss+1)*Eps_T(1,1,qq,ss+1,0) -SG%Jacob_sqf(qq,ss+2)*Eps_T(1,1,qq,ss+2,0)) &
                        +0.5*inv_du(2)*(SG%Jacob_sqf(qq+1,ss)*Eps_T(2,1,qq+1,ss,0) -SG%Jacob_sqf(qq-1,ss)*Eps_T(2,1,qq-1,ss,0)) &
                        +0.5*inv_du(3)*(SG%Jacob_sqf(qq,ss)*Eps_T(3,1,qq,ss,+1) -SG%Jacob_sqf(qq,ss)*Eps_T(3,1,qq,ss,-1));

               val= -(-2.0*Eps_T(1,1,qq,ss,0)*inv_du(1) + Kss_coeff/(SG%Jacob_sqf(qq,ss)))
               val=val*BC0_value ! multiply dfds @s=s0
               call VecSetValues(bVec,1,Gid,val,ADD_VALUES,ierr); CHKERRQ(ierr)
            endif
         enddo
      enddo
   endif

   call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
   call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)


   if(BC1_type==0) then ! Dirichlet condition (phi at outer boundary)
      ss=dim_a
      do ff=1,mzetamax
         do qq=1,dim_gq-1
            Lid=mapP2L(qq,ss,ff)
            if(Lid>0) then
               Gid=mapL2G(Lid)
               val=BC1_value
               call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
            endif
         enddo
      enddo
   elseif(BC1_type==1) then
      ss=dim_a
      inv_du(1)=-1.0*SG%inv_ds ! since delta_s is now negative value for outer boundary case
      do ff=1,mzetamax
         do qq=1,dim_gq-1
            Lid=mapP2L(qq,ss,ff)
            if(Lid>0) then
               Gid=mapL2G(Lid)
               Kss_coeff = 0.5*inv_du(1)*(-3*SG%Jacob_sqf(qq,ss)*Eps_T(1,1,qq,ss,0) +4*SG%Jacob_sqf(qq,ss-1)*Eps_T(1,1,qq,ss-1,0) -SG%Jacob_sqf(qq,ss-2)*Eps_T(1,1,qq,ss-2,0)) &
                        +0.5*inv_du(2)*(SG%Jacob_sqf(qq+1,ss)*Eps_T(2,1,qq+1,ss,0) -SG%Jacob_sqf(qq-1,ss)*Eps_T(2,1,qq-1,ss,0)) &
                        +0.5*inv_du(3)*(SG%Jacob_sqf(qq,ss)*Eps_T(3,1,qq,ss,+1) -SG%Jacob_sqf(qq,ss)*Eps_T(3,1,qq,ss,-1));

               val= -(-2.0*Eps_T(1,1,qq,ss,0)*inv_du(1) + Kss_coeff/(SG%Jacob_sqf(qq,ss)))
               val=val*BC1_value ! multiply dfds @s=s0
               call VecSetValues(bVec,1,Gid,val,ADD_VALUES,ierr); CHKERRQ(ierr)
            endif
         enddo
      enddo
   endif

   call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
   call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)


   if(BC0_type==1 .and. BC1_type==1) then
      ! This is purely Neumann conditions.. so the solution is not unique
      ! We need a reference point (Dirichlet condition)
      Lid=mapP2L(1,1,1) ! This is a reference point
      if(Lid>0) then
         Gid=mapL2G(Lid)
         val=0.0 ! Reference value
         call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      endif
   endif

   call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
   call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)


end subroutine Apply_BC_to_GK_Poisson_RHS_Coeffs

end module MGY_PETSc_Interface_Data