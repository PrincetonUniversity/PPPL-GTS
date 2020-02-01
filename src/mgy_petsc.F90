module MGY_PETSc_Interface_Data
   #include <petsc/finclude/petsc.h>
   use petscksp
   implicit none
   type PEST_Coord
      PetscInt gqid, aid, gfid
      PetscReal gq, a, gf
   end type PEST_Coord


   type(PEST_Coord),allocatable,dimension(:) :: mapG2P ! From Global node index to PEST node indices
   PetscInt, allocatable, dimension(:,:,:) :: mapP2G, mapP2L ! From PEST indices to Global or Local node index
   PetscInt, allocatable, dimension(:) :: mapL2G ! From Local indices to Global indices (this depends on the MPI process)
   PetscInt Lsize, Gsize, sum_Lsize_partd ! Size of Local nodes and Global nodes
   PetscInt, allocatable, dimension(:) :: recv_Lsize_partd, displs_Lsize_partd
   PetscInt, allocatable, dimension(:) :: recv_mapL2G_partd
   PetscInt Gid_start, Gid_end
   Mat Amat
   Vec bVec, xVec, ResVec
   KSP ksp

end module MGY_PETSc_Interface_Data

subroutine Initialize_MGY_PETSc_Interface_Data(dim_gq,dim_a,dim_gf)
   use MGY_PETSc_Interface_Data
   use precision
   use global_parameters,only:stdout,mype,numberpe
   use particle_decomp
   implicit none
   integer, intent(in):: dim_gq,dim_a,dim_gf
   integer :: gqid, aid, gfid, Gid, Lid, tnum, extra, i
   integer, allocatable :: d_nnz(:), o_nnz(:)
   PetscInt neq_All3D, neq_Pol2D
   integer, allocatable :: tarGid(:,:)
   PetscReal, allocatable :: tarVal(:,:)
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
   neq_All3D = Gsize

   ! Manual partitioning for each poloidal plane
   neq_Pol2D = (dim_gq-1)*dim_a
   Lsize = neq_Pol2D/npartdom
   extra = neq_Pol2D-Lsize*npartdom

   Gid_start = myrank_toroidal*neq_Pol2D + Lsize*myrank_partd + min(myrank_partd,extra)
   Gid_end   = myrank_toroidal*neq_Pol2D + Lsize*(myrank_partd+1) + min(myrank_partd+1,extra)-1
   Lsize = Gid_end-Gid_start+1


   allocate(mapL2G(Lsize))
   allocate(mapG2P(0:neq_All3D-1)) ! 1-dimensional array (starting index is 0 because of PETSc's C-style indexing)
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

   mapP2L=-100 ! Default value is negative (-100), so one can tell whether a specific PETS point is out of the its process (mapP2L<0)
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

   ! Calculate target Gid & value arrays which will be used to construct Amat
   allocate(tarGid(7,Lsize), tarVal(7,Lsize))
   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   do Lid=1,Lsize
      Gid=mapL2G(Lid)

      gqid=mapG2P(Gid)%gqid
      aid =mapG2P(Gid)%aid
      gfid=mapG2P(Gid)%gfid

      tnum=1
      tarGid(tnum,Lid)=mapP2G(gqid,aid,gfid);   tarVal(tnum,Lid)=6.0;  tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid-1,aid,gfid); tarVal(tnum,Lid)=-1.0; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid+1,aid,gfid); tarVal(tnum,Lid)=-1.0; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid,aid-1,gfid); tarVal(tnum,Lid)=-1.0; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid,aid+1,gfid); tarVal(tnum,Lid)=-1.0; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid,aid,gfid-1); tarVal(tnum,Lid)=-1.0; tnum=tnum+1;
      tarGid(tnum,Lid)=mapP2G(gqid,aid,gfid+1); tarVal(tnum,Lid)=-1.0; tnum=tnum+1;
   enddo

   ! To calculate numbers of non-zero elements for diagonal(d_nnz) & off-diagonal(o_nnz) blocks
   allocate(d_nnz(Lsize), o_nnz(Lsize))
   d_nnz=0; o_nnz=0;
   do Lid=1,Lsize
      do tnum=1,7
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
   call MatSetSizes(Amat,Lsize,Lsize,neq_All3D,neq_All3D,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)
   ! Preallocate the memory using o_nnz and d_nnz
   call MatMPIAIJSetPreallocation(Amat,7,d_nnz,7,o_nnz,ierr); CHKERRQ(ierr)
   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(Amat,1,Gid,7,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo


   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)


   call VecCreate(PETSC_COMM_WORLD,bVec,ierr); CHKERRQ(ierr)
   call VecSetSizes(bVec,Lsize,neq_All3D,ierr); CHKERRQ(ierr)
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

   ! call VecSet(bVec,5.0,ierr); CHKERRQ(ierr)
   ! call VecGetArrayF90(bVec,ptr_bVec,ierr)
   ! ptr_bVec(1) = 123.0
   ! call VecRestoreArrayF90(bVec,ptr_bVec,ierr)
   ! call VecView(bVec,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   write(stdout,"(3(A,I8))") "mype=",mype,", Gid_start=",Gid_start,", Gid_end=",Gid_end
   call flush(stdout)

   write(stdout,"(2(A,I8),1x,I8,A,I8)") "Gsize=",Gsize,", Lsize=",Lsize, 7*Lsize, ", sum(d_nnz+o_nnz)=",sum(d_nnz+o_nnz)
   call flush(stdout)
   ! call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

   deallocate(tarGid,tarVal)
   deallocate(d_nnz,o_nnz)


end subroutine Initialize_MGY_PETSc_Interface_Data

subroutine Finalize_MGY_PETSc_Interface_Data
   use MGY_PETSc_Interface_Data
   implicit none

   deallocate(mapG2P,mapP2G,mapP2L,mapL2G)
   deallocate(recv_Lsize_partd,displs_Lsize_partd)
   deallocate(recv_mapL2G_partd)

end subroutine Finalize_MGY_PETSc_Interface_Data


subroutine Construct_Poisson3D_Solver_using_DiagnosticGrid
   use MGY_PETSc_Interface_Data
   use Diagnostic_Grid
   use particle_decomp
   use global_parameters,only:mype,stdout,zetamax,zetamin,mzetamax
   implicit none

   PetscInt Lid,Gid,gqid,aid,gfid,tnum
   PetscInt, allocatable :: tarGid(:,:)
   PetscReal Coeff_gq, Coeff_a, Coeff_gf
   PetscReal, allocatable :: tarVal(:,:)
   real tt0, tt1
   PetscErrorCode ierr
   real, external:: dsda

   call Initialize_MGY_PETSc_Interface_Data(DG%dim_gq,DG%dim_a,ntoroidal)

   allocate(tarGid(7,Lsize), tarVal(7,Lsize))

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

   ! Apply Neumann Condition at inner and outer boundary
   do gfid=1,ntoroidal
      do gqid=1,DG%dim_gq-1
         ! Neumann condition (Er=0 at inner, outer boundary)
         Lid=mapP2L(gqid,1,gfid)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=-1.0 !phi(aid=1)
            tarVal(5,Lid)=+1.0 !phi(aid=2)
         endif
         Lid=mapP2L(gqid,DG%dim_a,gfid)
         if(Lid>0) then
            tarVal(:,Lid)=0
            tarVal(1,Lid)=-1.0 !phi(aid=dim_a)
            tarVal(4,Lid)=+1.0 !phi(aid=dim_a-1)
         endif
      enddo
   enddo
   ! Set Reference potential value 0 at (gqid,aid,gfid)=(1,1,1)
   if(mapP2L(1,2,1)>0) then
      Lid=mapP2L(1,2,1)
      tarVal(:,Lid)=0
      tarVal(1,Lid)=1
   endif


   ! if(mype==0) then
   !    do Lid=1,Lsize
   !       Gid=mapL2G(Lid)
   !       gqid=mapG2P(Gid)%gqid
   !       aid=mapG2P(Gid)%aid
   !       gfid=mapG2P(Gid)%gfid
   !       write(stdout,"(2(A,I5),A,3(I5,1x),A,7(ES10.2,1x))") "Lid=",Lid,",Gid=",Gid,",(gq,a,gf)id=(",gqid,aid,gfid, "), tarVal=", tarVal(:,Lid)
   !    enddo
   !    call flush(stdout)
   ! endif

   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      call MatSetValues(Amat,1,Gid,7,tarGid(:,Lid),tarVal(:,Lid),INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   ! Assembly Matrix
   tt0=MPI_WTIME()
   call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "After MatSetValues using DG, MatAssembly takes ",tt1-tt0," secs"
   call flush(stdout)


   ! Construct KSP
   tt0=MPI_WTIME()
   call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp,Amat,Amat,ierr); CHKERRQ(ierr)
   call KSPSetTolerances(ksp,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPCreate takes ",tt1-tt0," secs"
   call flush(stdout)

   deallocate(tarGid, tarVal)

end subroutine Construct_Poisson3D_Solver_using_DiagnosticGrid

subroutine Solve_Poisson3D_using_DiagnosticGrid()
   use precision, only:wp
   use global_parameters, only:stdout,mype,ifelectron,mzetamax
   use MGY_PETSc_Interface_Data
   use Diagnostic_Grid
   use particle_decomp
   implicit none

   PetscScalar, pointer :: ptr_xVec(:), ptr_bVec(:), ptr_ResVec(:)
   PetscInt Lid,Gid,its, gqid,gfid
   PetscReal val, res_norm2, src_norm2
   PetscErrorCode ierr
   PetscReal source(DG%dim_gq,DG%dim_a)
   PetscReal, allocatable :: xVec_mype(:),xVec_partd(:),bVec_mype(:),bVec_partd(:), ResVec_mype(:),ResVec_partd(:)
   real(wp) tt0,tt1

   source = DG%iScaData(:,:,1,1)-(DG%eScaData(:,:,1,1)+int(ifelectron/2)*DG%dne_adiab)

   ! Assign and Assemble bVec (source term)
   do Lid=1,Lsize
      Gid=mapL2G(Lid)
      val=source(mapG2P(Gid)%gqid,mapG2P(Gid)%aid)
      call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
   enddo

   ! Set Zero bValues at both inner&outer boundary fo Neumann Condition
   do gfid=1,mzetamax
      do gqid=1,DG%dim_gq-1
         ! Neumann condition (Er=0 at inner, outer boundary)
         Lid=mapP2L(gqid,1,gfid)
         if(Lid>0) then
            Gid=mapL2G(Lid)
            val=0
            call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
         endif
         Lid=mapP2L(gqid,DG%dim_a,gfid)
         if(Lid>0) then
            Gid=mapL2G(Lid)
            val=0
            call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
         endif
      enddo
   enddo

   ! Set Zero bValue at the reference value
   if(mapP2L(1,2,1)>0) then
      Lid=mapP2L(1,2,1)
      Gid=mapL2G(Lid)
      val=0
      call VecSetValues(bVec,1,Gid,val,INSERT_VALUES,ierr); CHKERRQ(ierr)
   endif




   call VecAssemblyBegin(bVec,ierr); CHKERRQ(ierr)
   call VecAssemblyEnd(bVec,ierr); CHKERRQ(ierr)

   ! Solve the linear system
   tt0=MPI_WTIME()
   call KSPSolve(ksp,bVec,xVec,ierr); CHKERRQ(ierr)
   tt1=MPI_WTIME()
   if(mype==0) write(stdout,*) "KSPSolve takes ",tt1-tt0," secs"
   call flush(stdout)

   ! Check the residual error and number of the iterations
   call MatMult(Amat,xVec,ResVec,ierr); CHKERRQ(ierr)
   call VecAXPY(ResVec,-1.0,bVec,ierr); CHKERRQ(ierr)
   call VecNorm(ResVec,NORM_2,res_norm2,ierr); CHKERRQ(ierr)
   call VecNorm(bVec,NORM_2,src_norm2,ierr); CHKERRQ(ierr)
   call KSPGetIterationNumber(ksp,its,ierr); CHKERRQ(ierr)
   if(mype==0) write(stdout,"(A,ES16.8,A,4I)") "Norm of error=", res_norm2/src_norm2, ", iterations=",its
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
   PetscInt :: gid, fGidx, MGidx_start, MGidx_end, VGidx_start, VGidx_end, nx=30,ny=30, neq, idxm, its
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
   call KSPGetIterationNumber(ksp,its,ierr); CHKERRQ(ierr)

   if(mype==0) write(stdout,"(A,ES16.8,A,4I)") "Norm of error=", x_norm2, ", iterations=",its

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



! subroutine test_Poisson_3D()
!    use MGY_PETSc_Interface_Data
!    use Diagnostic_Grid
!    use global_parameters
!    implicit none

!    Vec x,b,u
!    Mat Amat
!    KSP ksp
!    PetscErrorCode ierr
!    PetscInt :: gid, fGidx, MGidx_start, MGidx_end, VGidx_start, VGidx_end, nx=200*16,ny=200*16, neq, idxm, its
!    PetscInt :: mype, numberpe

!    call Initialize_MGY_PETSc_Interface_Data(DG%dim_gq,DG%dim_a,mzetamax)

!    neq = Gsize

!    call PetscInitialize('./PETSc_input',ierr); CHKERRQ(ierr)
!    call mpi_comm_size(mpi_comm_world,numberpe,ierr)
!    call mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)

!    call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHKERRQ(ierr)
!    call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,neq,neq,ierr); CHKERRQ(ierr)
!    call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)
!    call MatMPIAIJSetPreallocation(Amat,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,ierr); CHKERRQ(ierr)






!    call Finalize_MGY_PETSc_Interface_Data


! end subroutine test_Poisson_3D
