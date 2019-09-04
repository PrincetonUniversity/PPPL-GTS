#include <stdio.h>
#include "gts_git_diff.h"
#include "gts_module_files.h"
#define Print_Git_Info print_git_info_

int Print_Git_Info(){

#ifdef GIT_VERSION
   printf("============ GTS VERSION INFORMATION ============\n");
   printf("\n");
   printf("       VERSION: %s\n", GIT_VERSION);
   printf("        BRANCH: %s\n", GIT_BRANCH);
   printf("   COMMIT  MSG: %s\n", GIT_COM_MSG);
   printf("   COMMIT  SHA: %s\n", GIT_COM_SHA);
   printf("   COMMIT DATE: %s\n", GIT_DATE);
   printf("  WORKING TREE: %s\n", GIT_WTREE);
   printf("\n");
   printf("    BUILD DATE: %s\n", BUILD_DATE);
   printf("    BUILD  DIR: %s\n", BUILD_DIR);
   printf("\n");
   printf("=================================================\n");
   printf("\n");
   printf("If \"WORKING TREE\" is not clean, check \"GTS_GIT_DIFF.out\" for details\n");
   printf("\n");
   printf("\n");
   printf("%s\n",gts_module_files_out);
   fflush(stdout);

   FILE *fid2 = fopen("GTS_VERSION_INFO.out","wt");
   fprintf(fid2,"============ GTS VERSION INFORMATION ============\n");
   fprintf(fid2,"\n");
   fprintf(fid2,"       VERSION: %s\n", GIT_VERSION);
   fprintf(fid2,"        BRANCH: %s\n", GIT_BRANCH);
   fprintf(fid2,"   COMMIT  MSG: %s\n", GIT_COM_MSG);
   fprintf(fid2,"   COMMIT  SHA: %s\n", GIT_COM_SHA);
   fprintf(fid2,"   COMMIT DATE: %s\n", GIT_DATE);
   fprintf(fid2,"  WORKING TREE: %s\n", GIT_WTREE);
   fprintf(fid2,"\n");
   fprintf(fid2,"    BUILD DATE: %s\n", BUILD_DATE);
   fprintf(fid2,"    BUILD  DIR: %s\n", BUILD_DIR);
   fprintf(fid2,"\n");
   fprintf(fid2,"=================================================\n");
   fprintf(fid2,"\n");
   fprintf(fid2,"If \"WORKING TREE\" is not clean, check \"GTS_GIT_DIFF.out\" for details\n");
   fprintf(fid2,"\n");
   fprintf(fid2,"\n");
   fprintf(fid2,"%s\n",gts_module_files_out);
   fclose(fid2);
#endif

   FILE *fid= fopen("GTS_GIT_DIFF.out","w+");
   fprintf(fid,"%s\n",gts_git_diff_out);
   fclose(fid);

   return(0);
}
