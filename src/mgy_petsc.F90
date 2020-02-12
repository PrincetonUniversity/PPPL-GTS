module MGY_PETSc_Interface_Data
   #include <petsc/finclude/petsc.h>
   use petscksp
   use precision, only:wp
   implicit none
   type PEST_Indices
      PetscInt gqid, aid, gfid
   end type PEST_Indices

   type(PEST_Indices),allocatable,dimension(:) :: mapG2P ! From Global node index to PEST node indices
   PetscInt, allocatable, dimension(:,:,:) :: mapP2G, mapP2L ! From PEST indices to Global or Local node index
   PetscInt, allocatable, dimension(:) :: mapL2G ! From Local indices to Global indices (this depends on the MPI process)
   PetscInt Lsize, Gsize, sum_Lsize_partd ! Size of Local nodes and Global nodes
   PetscInt, allocatable, dimension(:) :: recv_Lsize_partd, displs_Lsize_partd
   PetscInt, allocatable, dimension(:) :: recv_mapL2G_partd
   PetscInt Gid_start, Gid_end
   Mat Amat
   Vec bVec, xVec, ResVec
   KSP ksp

   PetscInt, allocatable, dimension(:,:) :: tarGid
   PetscReal, allocatable, dimension(:,:) :: tarVal

   ! Following Tensors will be used to determine Amat components (Poisson equation)
   ! Note that they are composed with the contravariant components
   ! gg_T:  contravariant Metrics; g^(ij)=e^i*e^j, where e is covariant vector on the PEST coordinate
   ! polar_T: Polarization shielding tensor. (n_s*m_s/B^2)*(I_tensor-b*b), where b is unit vector of the magnetic field
   ! eps_T:
   ! CT:

   real(wp), allocatable, dimension(:,:,:,:) :: Metric_T, Polar_T ! On the poloidal plane (PEST structured coordinate)
   real(wp), allocatable, dimension(:,:,:,:,:) :: Eps_T, CT ! Including neighboring left and right poloidal planes


   ! Note that indices are (-1,0,1) for each direction. This is very convenient to describe neighboring points.
   ! Vec Metric_T(-1:1,-1:1), Polar_T(-1:1,-1:1), Eps_T(-1:1,-1:1), CT(-1:1,-1:1)

end module MGY_PETSc_Interface_Data

subroutine Initialize_MGY_PETSc_Interface_Data(dim_gq,dim_a,dim_gf)
   use MGY_PETSc_Interface_Data
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
   ! call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

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
               if(mapP2G(gqid,aid,gfid)>0) then
                  write(200205,"(5(I8))") mapP2L(gqid,aid,gfid), mapP2G(gqid,aid,gfid), gqid, aid, gfid
               endif
            enddo
         enddo
      enddo
      call flush(200205)
      close(200205)

   endif


end subroutine Initialize_MGY_PETSc_Interface_Data

subroutine Finalize_MGY_PETSc_Interface_Data
   use MGY_PETSc_Interface_Data
   implicit none

   deallocate(mapG2P,mapP2G,mapP2L,mapL2G)
   deallocate(recv_Lsize_partd,displs_Lsize_partd)
   deallocate(recv_mapL2G_partd)

   if(allocated(Metric_T)) deallocate(Metric_T, Polar_T, Eps_T)
   if(allocated(CT)) deallocate(CT)

end subroutine Finalize_MGY_PETSc_Interface_Data


subroutine Construct_SimplePoisson3D_Solver_using_DiagnosticGrid
   use MGY_PETSc_Interface_Data
   use Diagnostic_Grid
   use particle_decomp
   use global_parameters,only:mype,stdout,zetamax,zetamin,mzetamax
   implicit none

   PetscInt Lid,Gid,gqid,aid,gfid,tnum
   PetscReal Coeff_gq, Coeff_a, Coeff_gf
   real tt0, tt1
   PetscErrorCode ierr
   real, external:: dsda
   integer :: max_nnz=7
   integer, allocatable :: d_nnz(:),o_nnz(:)

   call Initialize_MGY_PETSc_Interface_Data(DG%dim_gq,DG%dim_a,ntoroidal)

   allocate(tarGid(max_nnz,Lsize), tarVal(max_nnz,Lsize))

   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   do Lid=1,Lsize
      Gid=mapL2G(Lid)

      gqid=mapG2P(Gid)%gqid
      aid =mapG2P(Gid)%aid
      gfid=mapG2P(Gid)%gfid

      Coeff_gq = DG%g_qq_cntv(gqid,aid)*DG%inv_dgq**2.0
      Coeff_a  = DG%g_aa_cntv(gqid,aid)*(dsda(DG%a_1D(aid))*DG%inv_ds)**2.0
      Coeff_gf = DG%g_ff_cntv(gqid,aid)/((zetamax-zetamin)**2.0)

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
      do gqid=1,DG%dim_gq-1
         ! Neumann condition (Er=0 at inner, outer boundary)
         Lid=mapP2L(gqid,2,gfid)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=-1.0 !phi(aid=2)
            tarVal(4,Lid)=+1.0 !phi(aid=1)
         endif
         Lid=mapP2L(gqid,DG%dim_a-1,gfid)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=-1.0 !phi(aid=dim_a-1)
            tarVal(5,Lid)=+1.0 !phi(aid=dim_a)
         endif

         ! Neumann condition (E_theta=0 at inner, outer boundary)
         Lid=mapP2L(gqid,DG%dim_a,gfid)
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

   ! Construct Matrix (Amat)
   call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHKERRQ(ierr)
   call MatSetSizes(Amat,Lsize,Lsize,Gsize,Gsize,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)
   ! Preallocate the memory using o_nnz and d_nnz
   call MatMPIAIJSetPreallocation(Amat,max_nnz,d_nnz,max_nnz,o_nnz,ierr); CHKERRQ(ierr)
   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(Amat,1,Gid,max_nnz,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo


   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

   call VecCreate(PETSC_COMM_WORLD,bVec,ierr); CHKERRQ(ierr)
   call VecSetSizes(bVec,Lsize,Gsize,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(bVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,xVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,ResVec,ierr); CHKERRQ(ierr)


   ! Construct KSP
   tt0=MPI_WTIME()
   call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp,Amat,Amat,ierr); CHKERRQ(ierr)
   call KSPSetTolerances(ksp,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
   call flush(stdout)

   deallocate(d_nnz, o_nnz)

end subroutine Construct_SimplePoisson3D_Solver_using_DiagnosticGrid


subroutine Solve_Poisson3D_using_DiagnosticGrid()
   use precision, only:wp
   use global_parameters, only:stdout,mype,ifelectron,mzetamax
   use MGY_PETSc_Interface_Data
   use Diagnostic_Grid
   use particle_decomp
   use field_array, only: test_src
   implicit none

   PetscScalar, pointer :: ptr_xVec(:), ptr_bVec(:), ptr_ResVec(:)
   PetscInt Lid,Gid,num_its, gqid,gfid,aid
   PetscInt ss,qq,ff
   PetscReal val, res_norm2, src_norm2, RHS_leading_coeff, Kss_coeff,inv_du(3)
   PetscErrorCode ierr
   PetscReal source(DG%dim_gq,DG%dim_a)
   PetscReal, allocatable :: xVec_mype(:),xVec_partd(:),bVec_mype(:),bVec_partd(:), ResVec_mype(:),ResVec_partd(:)
   real(wp) tt0,tt1
   real(wp), external :: density

   ! source = DG%iScaData(:,:,1,1)-(DG%eScaData(:,:,1,1)+int(ifelectron/2)*DG%dne_adiab)
   source=0.0
   ! do aid=20,30
   !    do gqid=1,DG%dim_gq
   !       source(gqid,aid)=sin(6*(DG%gq_1D(gqid))-4*(DG%zetamin))
   !    enddo
   ! enddo
   call Interpolate_OriGridData_into_DG(test_src,source)

   RHS_leading_coeff=1.602176634e-19/8.8541878128e-12  ! element_charge/permittivity (e/eps0)

   source=RHS_leading_coeff*source
   do aid=1,DG%dim_a
      source(:,aid) = 1e6*density(DG%a_1D(aid),0.0)*source(:,aid)
   enddo

   ! Set Zero bValues at both inner&outer boundary fo Neumann Condition

   ! Assign and Assemble bVec (source term)
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      val=source(mapG2P(Gid)%gqid,mapG2P(Gid)%aid)
      call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   inv_du(1) = DG%inv_ds
   inv_du(2) = DG%inv_dgq
   inv_du(3) = 1.0/(DG%zetamax-DG%zetamin)


   ! Apply  Dirichlet Condition at outer boundary
   do ff=1,mzetamax
      do qq=1,DG%dim_gq-1
         ! Dirichlet condition (phi=0 at outer boundary)
         ss=DG%dim_a
         Lid=mapP2L(qq,ss,ff)
         if(Lid>0) then
            Gid=mapL2G(Lid)
            val=0
            call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
         endif

         ! Neumann condition (Er value at inner boundary)
         ss=1
         Lid=mapP2L(qq,ss,ff)
         if(Lid>0) then
            Gid=mapL2G(Lid)
            Kss_coeff =    inv_du(1)*(-3*DG%Jacob_sqf(qq,ss)*Eps_T(1,1,qq,ss,0) +4*DG%Jacob_sqf(qq,ss+1)*Eps_T(1,1,qq,ss+1,0) -DG%Jacob_sqf(qq,ss+2)*Eps_T(1,1,qq,ss+2,0)) &
                      +0.5*inv_du(2)*(DG%Jacob_sqf(qq+1,ss)*Eps_T(2,1,qq+1,ss,0) -DG%Jacob_sqf(qq-1,ss)*Eps_T(2,1,qq-1,ss,0)) &
                      +0.5*inv_du(3)*(DG%Jacob_sqf(qq,ss)*Eps_T(3,1,qq,ss,+1) -DG%Jacob_sqf(qq,ss)*Eps_T(3,1,qq,ss,-1));

            val= -(-2.0*Eps_T(1,1,qq,ss,00)*inv_du(1) + Kss_coeff/(DG%Jacob_sqf(qq,ss)))
            val=val*0.0 ! multiply dfds @s=s0
            call VecSetValues(bVec,1,Gid,val,ADD_VALUES,ierr); CHKERRQ(ierr)
         endif
      enddo
   enddo





   call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
   call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)

   ! Solve the linear system
   tt0=MPI_WTIME()
   call KSPSolve(ksp,bVec,xVec,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPSolve takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD, ierr);CHKERRQ(ierr);

   ! Check the residual error and number of the iterations
   call MatMult(Amat,xVec,ResVec,ierr); CHKERRQ(ierr)
   call VecAXPY(ResVec,-1.0,bVec,ierr); CHKERRQ(ierr)
   call VecNorm(ResVec,NORM_2,res_norm2,ierr); CHKERRQ(ierr)
   call VecNorm(bVec,NORM_2,src_norm2,ierr); CHKERRQ(ierr)
   call KSPGetIterationNumber(ksp,num_its,ierr); CHKERRQ(ierr)
   if(mype==0) write(stdout,"(A,ES16.8,A,4I)") "Norm of error=", res_norm2/src_norm2, ", iterations=",num_its
   call flush(stdout)

   ! Copy xVec (solution) to the DG%test_phi vector
   allocate(xVec_mype(Lsize), xVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(xVec,ptr_xVec,ierr)
   xVec_mype(:)=ptr_xVec(1:Lsize)
   call VecRestoreArrayReadF90(xVec,ptr_xVec,ierr)
   call MPI_ALLGATHERV(xVec_mype,Lsize,mpi_Rsize,xVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      DG%test_phi_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = xVec_partd(Lid)
   enddo
   DG%test_phi_2D(DG%dim_gq,:)=DG%test_phi_2D(1,:)
   deallocate(xVec_mype, xVec_partd)


   ! Copy bVec (source) to the DG%test_phi vector
   allocate(bVec_mype(Lsize), bVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(bVec,ptr_bVec,ierr)
   bVec_mype(:)=ptr_bVec(1:Lsize)
   call VecRestoreArrayReadF90(bVec,ptr_bVec,ierr)
   call MPI_ALLGATHERV(bVec_mype,Lsize,mpi_Rsize,bVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      DG%test_src_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = bVec_partd(Lid)
   enddo
   DG%test_src_2D(DG%dim_gq,:)=DG%test_src_2D(1,:)
   deallocate(bVec_mype, bVec_partd)

   ! Copy ResVec (residual=A*x-b) to the DG%test_phi vector
   allocate(ResVec_mype(Lsize), ResVec_partd(sum_Lsize_partd))
   call VecGetArrayReadF90(ResVec,ptr_ResVec,ierr)
   ResVec_mype(:)=ptr_ResVec(1:Lsize)
   call VecRestoreArrayReadF90(ResVec,ptr_ResVec,ierr)
   call MPI_ALLGATHERV(ResVec_mype,Lsize,mpi_Rsize,ResVec_partd,recv_Lsize_partd,displs_Lsize_partd,mpi_Rsize,partd_comm,ierr)
   do Lid=1,sum_Lsize_partd
      Gid=recv_mapL2G_partd(Lid)
      DG%test_res_2D(mapG2P(Gid)%gqid,mapG2P(Gid)%aid) = ResVec_partd(Lid)
   enddo
   DG%test_res_2D(DG%dim_gq,:)=DG%test_res_2D(1,:)
   deallocate(ResVec_mype, ResVec_partd)


   ! To get full 3D data for diagnostic purpose
   call MPI_ALLGATHER(DG%test_phi_2D,DG%dim_gq*DG%dim_a,mpi_Rsize,DG%test_phi_3D,DG%dim_gq*DG%dim_a,mpi_Rsize,toroidal_comm,ierr)
   call MPI_ALLGATHER(DG%test_src_2D,DG%dim_gq*DG%dim_a,mpi_Rsize,DG%test_src_3D,DG%dim_gq*DG%dim_a,mpi_Rsize,toroidal_comm,ierr)
   call MPI_ALLGATHER(DG%test_res_2D,DG%dim_gq*DG%dim_a,mpi_Rsize,DG%test_res_3D,DG%dim_gq*DG%dim_a,mpi_Rsize,toroidal_comm,ierr)

end subroutine Solve_Poisson3D_using_DiagnosticGrid


subroutine test_Laplacian_2D()
   #include <petsc/finclude/petsc.h>
   use petscksp
   use global_parameters,only:stdout,mype,numberpe
   implicit none

   Vec x,b,u
   Mat Amat
   KSP ksp
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

   call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHKERRQ(ierr)
   call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,neq,neq,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)
   ! Note that the following d_nz=5, o_nz=5 are not optimized values.
   call MatMPIAIJSetPreallocation(Amat,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,ierr); CHKERRQ(ierr)
   call MatGetOwnershipRange(Amat,MGidx_start,MGidx_end,ierr); CHKERRQ(ierr)

   MGidx_end=MGidx_end-1 !Note that this is an global index of the last local row more proper for Fortran do loop

   call PetscLogStageRegister("Assembly", stage,ierr); CHKERRQ(ierr)
   call PetscLogStagePush(stage,ierr); CHKERRQ(ierr)


      ! do gid=0,neq-1
      !    idxm=gid-1; val=-1.0
      !    if(idxm>=0) call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      !    idxm=gid+1; val=-1.0
      !    if(idxm<=neq-1) call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      !    idxm=gid-ny; val=-1.0
      !    if(idxm>=0) call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      !    idxm=gid+ny; val=-1.0
      !    if(idxm<=neq-1) call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)

      !    idxm=gid; val=4.0
      !    call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr)
      ! enddo

   do gid=MGidx_start,MGidx_end
      idxm=gid-1; val=-1.0
      if(idxm>=0) call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      idxm=gid+1; val=-1.0
      if(idxm<=neq-1) call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      idxm=gid-ny; val=-1.0
      if(idxm>=0) call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
      idxm=gid+ny; val=-1.0
      if(idxm<=neq-1) call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr); CHKERRQ(ierr)

      idxm=gid; val=4.0
      call MatSetValues(Amat,1,gid,1,idxm,val,INSERT_VALUES,ierr)
   enddo
   tt0=MPI_WTIME()
   call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)
   call PetscLogStagePop(ierr); CHKERRQ(ierr)

   ! call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

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
   call MatMult(Amat,u,b,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatMult takes ",tt1-tt0," secs"
   call flush(stdout)


   tt0=MPI_WTIME()
   call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp,Amat,Amat,ierr); CHKERRQ(ierr)
   call KSPSetTolerances(ksp,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
   call flush(stdout)

   tt0=MPI_WTIME()
   ! Solve the linear system
   call KSPSolve(ksp,b,x,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPSolve takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   ! Check solution and error
   call VecAXPY(x,-1.0,u,ierr); CHKERRQ(ierr)
   call VecNorm(x,NORM_2,x_norm2,ierr); CHKERRQ(ierr)
   call KSPGetIterationNumber(ksp,num_its,ierr); CHKERRQ(ierr)

   if(mype==0) write(stdout,"(A,ES16.8,A,4I)") "Norm of error=", x_norm2, ", iterations=",num_its

   call MPI_Barrier(PETSC_COMM_WORLD,ierr)
   write(stdout,"(A,I8,A,(I8,A,I8),A,(I8,A,I8),A)")"mype=",mype,", idx_start_(M,V)=(",MGidx_start,",",VGidx_start,"), idx_end_(M,V)=(",MGidx_end,",",VGidx_end,")"
   call flush(stdout)

   call KSPDestroy(ksp,ierr); CHKERRQ(ierr)
   call VecDestroy(u,ierr); CHKERRQ(ierr)
   call VecDestroy(x,ierr); CHKERRQ(ierr)
   call VecDestroy(b,ierr); CHKERRQ(ierr)
   call MatDestroy(Amat,ierr); CHKERRQ(ierr)

   ! call PetscFinalize(ierr)

   return
end subroutine test_Laplacian_2D


subroutine Construct_GK_Poisson3D_Solver_using_DiagnosticGrid
   use MGY_PETSc_Interface_Data
   use Diagnostic_Grid
   use particle_decomp
   use global_parameters,only:mype,stdout,zetamax,zetamin,mzetamax
   implicit none

   PetscInt Lid,Gid,qq,ss,ff,tnum
   PetscReal off_val
   real tt0, tt1
   PetscReal inv_Jacob
   PetscErrorCode ierr
   real, external:: dsda
   integer :: max_nnz=19 ! maximum number of non-zero values per each row
   integer, allocatable :: o_nnz(:), d_nnz(:)
   integer :: dim_gq, dim_a


   dim_a =DG%dim_a
   dim_gq=DG%dim_gq

   call Initialize_MGY_PETSc_Interface_Data(dim_gq,dim_a,ntoroidal)

   ! Tensors on the poloidal plane, and Coeff Tensor including neighboring left and right poloidal plane
   ! Note that there are ghost cells for easier indexing
   allocate(Metric_T(3,3,0:dim_gq,0:dim_a+1), Polar_T(3,3,0:dim_gq,0:dim_a+1))
   allocate(Eps_T(3,3,0:dim_gq,0:dim_a+1,-1:1), CT(3,3,0:dim_gq,0:dim_a+1,-1:1))
   call Determine_GK_PoissonTensors_on_DiagnosticGrid()

   allocate(tarGid(max_nnz,Lsize), tarVal(max_nnz,Lsize))

   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   do Lid=1,Lsize
      Gid=mapL2G(Lid)

      qq = mapG2P(Gid)%gqid
      ss = mapG2P(Gid)%aid
      ff = mapG2P(Gid)%gfid

      inv_Jacob=1.0/DG%Jacob_sqf(qq,ss)


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

   ! Apply Neumann (+Dirichlet) Condition at inner and outer boundary
   do ff=1,ntoroidal
      do qq=1,dim_gq-1
         ! ! Neumann condition (Er=0 at inner, outer boundary)
         ! Lid=mapP2L(qq,2,ff)
         ! if(Lid>0) then
         !    tarVal(:,Lid)=0
         !    tarVal(1,Lid)=-1.0 !phi(ss=2)
         !    tarVal(4,Lid)=+1.0 !phi(ss=1)
         ! endif
         ! Lid=mapP2L(qq,dim_a-1,ff)
         ! if(Lid>0) then
         !    tarVal(:,Lid)=0
         !    tarVal(1,Lid)=-1.0 !phi(ss=dim_a-1)
         !    tarVal(5,Lid)=+1.0 !phi(ss=dim_a)
         ! endif

         ! ! Neumann condition (E_theta=0 at inner, outer boundary)
         ! Lid=mapP2L(qq,dim_a,ff)
         ! if(Lid>0) then
         !    tarVal(:,Lid)=0
         !    tarVal(1,Lid)=-1.0 !phi(ss=dim_a, qq)
         !    tarVal(3,Lid)=+1.0 !phi(ss=dim_a, qq+1)
         ! endif

         !Dirichlet condition (phi=0 at inner boundary)
         Lid=mapP2L(qq,1,ff)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=1.0 !phi(ss=1)
         endif

         !Dirichlet condition (phi=0 at outer boundary)
         Lid=mapP2L(qq,dim_a,ff)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=1.0 !phi(ss=dim_a, qq)
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


   ! Construct Matrix (Amat)
   call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHKERRQ(ierr)
   call MatSetSizes(Amat,Lsize,Lsize,Gsize,Gsize,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)
   ! Preallocate the memory using o_nnz and d_nnz
   call MatMPIAIJSetPreallocation(Amat,max_nnz,d_nnz,max_nnz,o_nnz,ierr); CHKERRQ(ierr)
   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(Amat,1,Gid,max_nnz,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)

   call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)


   call VecCreate(PETSC_COMM_WORLD,bVec,ierr); CHKERRQ(ierr)
   call VecSetSizes(bVec,Lsize,Gsize,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(bVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,xVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,ResVec,ierr); CHKERRQ(ierr)

   ! Construct KSP
   tt0=MPI_WTIME()
   call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp,Amat,Amat,ierr); CHKERRQ(ierr)
   call KSPSetTolerances(ksp,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
   call flush(stdout)

   deallocate(d_nnz, o_nnz)

end subroutine Construct_GK_Poisson3D_Solver_using_DiagnosticGrid

subroutine Construct_GK_Poisson3D_Solver_with_Extremal_Method
   use MGY_PETSc_Interface_Data
   use Diagnostic_Grid
   use particle_decomp
   use global_parameters,only:mype,stdout,zetamax,zetamin,mzetamax
   implicit none

   PetscInt Lid,Gid,qq,ss,ff,tnum
   real tt0, tt1
   PetscReal inv_Jacob
   PetscErrorCode ierr
   real, external:: dsda
   integer :: max_nnz=19 ! maximum number of non-zero values per each row
   integer, allocatable :: o_nnz(:), d_nnz(:)
   integer :: dim_gq, dim_a
   PetscReal AntiSym1, AntiSym2, Sym1 ! anti-symmetric and symmetric terms

   dim_a =DG%dim_a
   dim_gq=DG%dim_gq

   call Initialize_MGY_PETSc_Interface_Data(dim_gq,dim_a,ntoroidal)

   ! Tensors on the poloidal plane, and Coeff Tensor including neighboring left and right poloidal plane
   ! Note that there are ghost cells for easier indexing
   allocate(Metric_T(3,3,0:dim_gq,0:dim_a+1), Polar_T(3,3,0:dim_gq,0:dim_a+1))
   allocate(Eps_T(3,3,0:dim_gq,0:dim_a+1,-1:1), CT(3,3,0:dim_gq,0:dim_a+1,-1:1))
   call Determine_GK_PoissonTensors_on_DiagnosticGrid()

   allocate(tarGid(max_nnz,Lsize), tarVal(max_nnz,Lsize))

   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   do Lid=1,Lsize
      Gid=mapL2G(Lid)

      qq = mapG2P(Gid)%gqid
      ss = mapG2P(Gid)%aid
      ff = mapG2P(Gid)%gfid

      inv_Jacob=1.0/DG%Jacob_sqf(qq,ss)

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
   call Apply_BC_to_Amat_Coeffs ! Apply Boundary Condition to Amt Coeffs





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


   ! Construct Matrix (Amat)
   call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHKERRQ(ierr)
   call MatSetSizes(Amat,Lsize,Lsize,Gsize,Gsize,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)
   ! Preallocate the memory using o_nnz and d_nnz
   call MatMPIAIJSetPreallocation(Amat,max_nnz,d_nnz,max_nnz,o_nnz,ierr); CHKERRQ(ierr)
   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(Amat,1,Gid,max_nnz,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)

   ! call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)


   call VecCreate(PETSC_COMM_WORLD,bVec,ierr); CHKERRQ(ierr)
   call VecSetSizes(bVec,Lsize,Gsize,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(bVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,xVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,ResVec,ierr); CHKERRQ(ierr)

   ! Construct KSP
   tt0=MPI_WTIME()
   call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp,Amat,Amat,ierr); CHKERRQ(ierr)
   call KSPSetTolerances(ksp,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
   call flush(stdout)

   deallocate(d_nnz, o_nnz)

end subroutine Construct_GK_Poisson3D_Solver_with_Extremal_Method


subroutine Determine_GK_PoissonTensors_on_DiagnosticGrid()
   use precision, only:wp
   use global_parameters, only:stdout,mype,ifelectron,mzetamax,msm0,amasse,amassi
   use MGY_PETSc_Interface_Data
   use Diagnostic_Grid
   use particle_decomp
   implicit none

   integer sid,gqid,i,j
   real(wp) invEpsilon,invBsq,n0,tmp ! 1/(B^2)
   real(wp), external :: density
   real(wp) inv_du(3)
   integer istatus(MPI_STATUS_SIZE), icount, ierror


   ! Metric Tensor is based on the static equilibrium magnetic field (PEST coordinate)
   do sid=1,DG%dim_a
      do gqid=1,DG%dim_gq
         Metric_T(:,:,gqid,sid)=0.0
         Metric_T(1,1,gqid,sid)=DG%g_ss_cntv(gqid,sid)
         Metric_T(2,1,gqid,sid)=DG%g_sq_cntv(gqid,sid)
         Metric_T(1,2,gqid,sid)=DG%g_sq_cntv(gqid,sid)
         Metric_T(2,2,gqid,sid)=DG%g_qq_cntv(gqid,sid)
         Metric_T(3,3,gqid,sid)=DG%g_ff_cntv(gqid,sid)
      enddo
   enddo

   ! Fill ghost nodes of Metric_T
   Metric_T(:,:,DG%dim_gq,:)=Metric_T(:,:,1,:)
   Metric_T(:,:,0,:)        =Metric_T(:,:,DG%dim_gq-1,:)
   Metric_T(:,:,:,0)        =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   Metric_T(:,:,:,DG%dim_a+1)  =-100 ! Assign Negative initial values which will be ignored in MatSetValues

   invEpsilon=1.0/8.8541878128e-12

   do sid=1,DG%dim_a
      ! Assuming the static equilibrium density for ions and electrons
      n0=1e6*density(DG%a_1D(sid),0) ! unit (#/m^(3))
      do gqid=1,DG%dim_gq
         invBsq=1.0/(DG%BB(gqid,sid)**2.0) ! unit (T^2)
         ! write(stdout,"(3(A,ES16.8))") "n0(sid)=",n0,", B=",DG%BB(gqid,sid), ", Polar(sid)=", invEpsilon*1e-3*(n0*amasse + n0*amassi)*invBsq
         do j=1,3
            do i=1,3
               tmp=invEpsilon*1e-3*(n0*amasse + n0*amassi)*invBsq  ! Polarization contribution from ions (n0*(1e-3*amassi)) and electrons (n0*1e-3*amasse)
               Polar_T(i,j,gqid,sid)= tmp*(Metric_T(i,j,gqid,sid)- invBsq*(DG%B0_cntv3_sqf(i,gqid,sid)*DG%B0_cntv3_sqf(j,gqid,sid)))
            enddo
         enddo
      enddo
   enddo
   Polar_T(:,:,DG%dim_gq,:)=Polar_T(:,:,1,:)
   Polar_T(:,:,0,:)        =Polar_T(:,:,DG%dim_gq-1,:)
   Polar_T(:,:,:,0)        =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   Polar_T(:,:,:,DG%dim_a+1)  =-100 ! Assign Negative initial values which will be ignored in MatSetValues

   ! Construct Epsilon Tensor ( eps=eps_0*Metric_T + Sigma_over_s(ns*ms/B^2 * (Metric_T - B^i*B^j/B^2)) )
   Eps_T(:,:,:,:,0)=Metric_T+Polar_T

   ! Communicate Eps_T with neighboring left and right planes.
   ! So, Eps_T(:,:,:,:,tid) ==>  tid=-1(left plane), tid=0(my plane), tid=1(right plane)
   ! Eps_T(3,3,0:dim_gq,0:dim_a+1,-1:1)
   icount=3*3*(DG%dim_gq+1)*(DG%dim_a+2)
   call MPI_SENDRECV(Eps_T(:,:,:,:,0),icount,mpi_Rsize,left_pe,0,Eps_T(:,:,:,:,1),icount,mpi_Rsize,right_pe,0,toroidal_comm,istatus,ierror)
   call MPI_SENDRECV(Eps_T(:,:,:,:,0),icount,mpi_Rsize,right_pe,0,Eps_T(:,:,:,:,-1),icount,mpi_Rsize,left_pe,0,toroidal_comm,istatus,ierror)



   inv_du(1) = DG%inv_ds
   inv_du(2) = DG%inv_dgq
   inv_du(3) = 1.0/(DG%zetamax-DG%zetamin)

   ! Determine Coefficient Tensor (CT = J*Eps_T(i,j)/(du^i*du^j))
   do sid=1,DG%dim_a
      do gqid=1,DG%dim_gq
         do j=1,3
            do i=1,3
               CT(i,j,gqid,sid,:) = DG%Jacob_sqf(gqid,sid)*Eps_T(i,j,gqid,sid,:)*inv_du(i)*inv_du(j)
            enddo
         enddo
      enddo
   enddo
   CT(:,:,DG%dim_gq,:,:)  =CT(:,:,1,:,:)
   CT(:,:,0,:,:)          =CT(:,:,DG%dim_gq-1,:,:)
   CT(:,:,:,0,:)          =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   CT(:,:,:,DG%dim_a+1,:) =-100 ! Assign Negative initial values which will be ignored in MatSetValues


end subroutine Determine_GK_PoissonTensors_on_DiagnosticGrid



subroutine Apply_BC_to_Amat_Coeffs
   use MGY_PETSc_Interface_Data
   use Diagnostic_Grid
   use global_parameters, only:mzetamax,stdout,mype,stdout
   implicit none

   PetscInt Lid,Gid,qq,ss,ff,j
   integer :: dim_gq, dim_a
   PetscReal inv_Jacob, Kj_coeff(3), inv_du(3)

   dim_a =DG%dim_a
   dim_gq=DG%dim_gq

   inv_du(1) = DG%inv_ds
   inv_du(2) = DG%inv_dgq
   inv_du(3) = 1.0/(DG%zetamax-DG%zetamin)

   ! Apply Neumann (+Dirichlet) Condition at inner and outer boundary
   do ff=1,mzetamax
      do qq=1,dim_gq-1
         ! Neumann condition (abs(E)=0 at inner)
         ss=1
         Lid=mapP2L(qq,ss,ff)
         if(Lid>0) then

            inv_Jacob=1.0/DG%Jacob_sqf(qq,ss)

            do j=1,3
               Kj_coeff(j) =    inv_du(1)*(-3*DG%Jacob_sqf(qq,ss)*Eps_T(1,j,qq,ss,0) +4*DG%Jacob_sqf(qq,ss+1)*Eps_T(1,j,qq,ss+1,0) -DG%Jacob_sqf(qq,ss+2)*Eps_T(1,j,qq,ss+2,0)) &
                           +0.5*inv_du(2)*(DG%Jacob_sqf(qq+1,ss)*Eps_T(2,j,qq+1,ss,0) -DG%Jacob_sqf(qq-1,ss)*Eps_T(2,j,qq-1,ss,0)) &
                           +0.5*inv_du(3)*(DG%Jacob_sqf(qq,ss)*Eps_T(3,j,qq,ss,+1) -DG%Jacob_sqf(qq,ss)*Eps_T(3,j,qq,ss,-1));
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



         ss=dim_a
         Lid=mapP2L(qq,ss,ff)
         if(Lid>0) then
            tarGid(:, Lid) = -100
            tarVal(:, Lid) = 0.0

            tarGid(1, Lid) = mapP2G(qq,ss,ff)
            tarVal(1, Lid) = 1.0 !phi(ss=dim_a)
         endif

      enddo
   enddo

end subroutine Apply_BC_to_Amat_Coeffs
