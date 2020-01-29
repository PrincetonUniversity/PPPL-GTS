subroutine test_Laplacian_2D()
   #include <petsc/finclude/petsc.h>
   use petscksp
   use global_parameters,only:stdout
   implicit none

   Vec x,b,u
   Mat Amat
   KSP ksp
   PetscReal val,x_norm2
   PetscInt :: gid, fGidx, MGidx_start, MGidx_end, VGidx_start, VGidx_end, nx=200,ny=100, neq, idxm, its, mype, numberpe
   PetscErrorCode ierr

   real tt0,tt1

#if defined(PETSC_USE_LOG)
   PetscLogStage stage;
#endif


   neq=nx*ny


   call PetscInitialize('./PETSc_input',ierr)
   if(ierr) return

   call mpi_comm_size(mpi_comm_world,numberpe,ierr)
   call mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
   if(mype==0) write(stdout,*) "numberpe=",numberpe
   call flush(stdout)

   call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHKERRQ(ierr)
   call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,neq,neq,ierr); CHKERRQ(ierr)
   call MatSetFromOptions(Amat,ierr); CHKERRQ(ierr)
   call MatMPIAIJSetPreallocation(Amat,5,PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,ierr); CHKERRQ(ierr)
   call MatGetOwnershipRange(Amat,MGidx_start,MGidx_end,ierr); CHKERRQ(ierr)

   MGidx_end=MGidx_end-1 !Note that this is an global index of the last local row more proper for Fortran do loop

   call PetscLogStageRegister("Assembly", stage,ierr); CHKERRQ(ierr)
   call PetscLogStagePush(stage,ierr); CHKERRQ(ierr)

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

   call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
   call PetscLogStagePop(ierr); CHKERRQ(ierr)

   call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)


   call VecCreate(PETSC_COMM_WORLD,u,ierr); CHKERRQ(ierr)
   call VecSetSizes(u,PETSC_DECIDE,neq,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(u,ierr); CHKERRQ(ierr)
   call VecDuplicate(u,b,ierr); CHKERRQ(ierr)
   call VecDuplicate(b,x,ierr); CHKERRQ(ierr)

   call VecGetOwnershipRange(u,VGidx_start,VGidx_end,ierr); CHKERRQ(ierr)
   VGidx_end=VGidx_end-1 !Note that this is an global index of the last local row more proper for Fortran do loop
   write(stdout,"(5(A,1x,I6))")"mype=",mype,", MGidx_start=",MGidx_start,", MGidx_end=",MGidx_end,", VGidx_start=",VGidx_start, ", VGidx_end=",VGidx_end
   call flush(stdout)

   call VecSet(u,1.0,ierr); CHKERRQ(ierr)

   ! call VecView(u,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   call MatMult(Amat,u,b,ierr); CHKERRQ(ierr)

   call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHKERRQ(ierr)
   call KSPSetOperators(ksp,Amat,Amat,ierr); CHKERRQ(ierr)
   call KSPSetTolerances(ksp,1e-5,PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr); CHKERRQ(ierr)
   call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)

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
   call flush(stdout)

   call KSPDestroy(ksp,ierr); CHKERRQ(ierr)
   call VecDestroy(u,ierr); CHKERRQ(ierr)
   call VecDestroy(x,ierr); CHKERRQ(ierr)
   call VecDestroy(b,ierr); CHKERRQ(ierr)
   call MatDestroy(Amat,ierr); CHKERRQ(ierr)

   call PetscFinalize(ierr)

   return
end subroutine test_Laplacian_2D