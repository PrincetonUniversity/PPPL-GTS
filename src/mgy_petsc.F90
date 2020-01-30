module MGY_PETSc_Interface_Data
   #include <petsc/finclude/petsc.h>
   use petscksp
   implicit none
   type PEST_Coord
      PetscInt gqid, aid, gfid
      PetscReal gq, a, gf
   end type PEST_Coord


   type(PEST_Coord),allocatable,dimension(:) :: mapG2P ! From Global node index to PEST node indices
   PetscInt, allocatable, dimension(:,:,:) :: mapP2G ! From PEST indices to Global node index
   PetscInt, allocatable, dimension(:) :: mapL2G ! From Local indices to Global indices (this depends on the MPI process)
   PetscInt num_of_Lnodes, num_of_Gnodes
   PetscInt Gid_start, Gid_end
   Mat Amat
   Vec bVec, xVec

end module MGY_PETSc_Interface_Data

subroutine Initialize_MGY_PETSc_Interface_Data(dim_gq,dim_a,dim_gf)
   use MGY_PETSc_Interface_Data
   use precision
   use global_parameters,only:stdout,mype,numberpe
   implicit none
   integer, intent(in):: dim_gq,dim_a,dim_gf
   integer :: gqid, aid, gfid, Gid, Lid, tnum
   integer, allocatable :: d_nnz(:), o_nnz(:)
   PetscInt neq
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



   num_of_Gnodes = (dim_gq-1)*dim_a*dim_gf
   neq = num_of_Gnodes

   allocate(mapG2P(0:neq-1)) ! 1-dimensional array (starting index is 0 because of PETSc's C-style indexing)
   allocate(mapP2G(0:dim_gq,0:dim_a+1,0:dim_gf+1)) ! 3-dimensional array (Note that 0, and dim+1 are ghost nodes)

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

   ! Fill ghost nodes of mapP2G
   mapP2G(dim_gq,:,:)   =mapP2G(1,:,:)
   mapP2G(0,:,:)        =mapP2G(dim_gq-1,:,:)
   mapP2G(:,0,:)        =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   mapP2G(:,dim_a+1,:)  =-100 ! Assign Negative initial values which will be ignored in MatSetValues
   mapP2G(:,:,0)        =mapP2G(:,:,dim_gf)
   mapP2G(:,:,dim_gf+1) =mapP2G(:,:,1)


   call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHKERRQ(ierr)
   call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,neq,neq,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)

   ! Preallocate the memory for matrix
   ! TODO : The following line is not optimal because it uses 7 non-zeros for both diagonal and off-diagonal
   !        This is not true. 7 is the number of non-zeros per each row (so diagonal+off-diagonal)
   !        The following will be optimized by determining and using d_nnz&o_nnz parameters for diagonal&off-diagonal,respectively
   call MatMPIAIJSetPreallocation(Amat,7,PETSC_NULL_INTEGER,7,PETSC_NULL_INTEGER,ierr); CHKERRQ(ierr)


   call MatGetOwnershipRange(Amat,Gid_start,Gid_end,ierr); CHKERRQ(ierr)
   Gid_end=Gid_end-1 ! For Fortran preferred indexing
   num_of_Lnodes=Gid_end-Gid_start+1

   allocate(mapL2G(num_of_Lnodes))
   do Lid=1,num_of_Lnodes
      mapL2G(Lid)=Gid_start+(Lid-1)
   enddo


   allocate(d_nnz(num_of_Lnodes), o_nnz(num_of_Lnodes))
   allocate(tarGid(7,num_of_Lnodes), tarVal(7,num_of_Lnodes))

   tarGid=-100 ! Assign Negative initial values which will be ignored in MatSetValues
   do Lid=1,num_of_Lnodes
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
   d_nnz=0; o_nnz=0;
   do Lid=1,num_of_Lnodes
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

   ! Preallocate the memory using o_nnz and d_nnz
   ! call MatMPIAIJSetPreallocation(Amat,7,d_nnz,7,o_nnz,ierr); CHKERRQ(ierr)
   call MatMPIAIJSetPreallocation(Amat,7,PETSC_NULL_INTEGER,7,PETSC_NULL_INTEGER,ierr); CHKERRQ(ierr)

   ! Set Mat Values using the non-zero information calculated from above
   do Lid=1,num_of_Lnodes
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
   call VecSetSizes(bVec,PETSC_DECIDE,neq,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(bVec,ierr); CHKERRQ(ierr)
   call VecDuplicate(bVec,xVec,ierr); CHKERRQ(ierr)

   call VecSet(bVec,5.0,ierr); CHKERRQ(ierr)

   ! call VecGetArrayF90(bVec,ptr_bVec,ierr)
   ! ptr_bVec(1) = 123.0
   ! call VecRestoreArrayF90(bVec,ptr_bVec,ierr)
   ! call VecView(bVec,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)


   write(stdout,"(2(A,I8),1x,I8,A,I8)") "num_of_Gnodes=",num_of_Gnodes,", num_of_Lnodes=",num_of_Lnodes, 7*num_of_Lnodes, ", sum(d_nnz+o_nnz)=",sum(d_nnz+o_nnz)
   call flush(stdout)
   ! call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)

   deallocate(tarGid,tarVal)
   deallocate(d_nnz,o_nnz)


end subroutine Initialize_MGY_PETSc_Interface_Data

subroutine Finalize_MGY_PETSc_Interface_Data
   use MGY_PETSc_Interface_Data
   implicit none

   deallocate(mapG2P,mapP2G,mapL2G)

end subroutine Finalize_MGY_PETSc_Interface_Data


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

   if(mype==0) write(stdout,*) "Norm of error=", x_norm2, ", iterations=",its

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

!    neq = num_of_Gnodes

!    call PetscInitialize('./PETSc_input',ierr); CHKERRQ(ierr)
!    call mpi_comm_size(mpi_comm_world,numberpe,ierr)
!    call mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)

!    call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHKERRQ(ierr)
!    call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,neq,neq,ierr); CHKERRQ(ierr)
!    call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)
!    call MatMPIAIJSetPreallocation(Amat,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,ierr); CHKERRQ(ierr)






!    call Finalize_MGY_PETSc_Interface_Data


! end subroutine test_Poisson_3D
