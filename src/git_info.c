#include <stdio.h>
#include "gts_git_diff.h"
#define Print_Git_Info print_git_info_

int Print_Git_Info(){


#ifdef GIT_VERSION
   printf("============ GTS VERSION INFORMATION ============\n");
   printf("\n");
   printf("       VERSION: %s\n", GIT_VERSION);
   printf("        BRANCH: %s\n", GIT_BRANCH);
   printf("    COMMIT MSG: %s\n", GIT_COM_MSG);
   printf("    COMMIT SHA: %s\n", GIT_COM_SHA);
   printf("      GIT DATE: %s\n", GIT_DATE);
   printf("    BUILD DATE: %s\n", BUILD_DATE);
   printf("  WORKING TREE: %s\n", GIT_WTREE);
   printf("\n");
   printf("=================================================\n");
   printf("\n");   
   printf("If \"WORKING TREE\" is not clean, check \"GTS_GIT_DIFF.out\" for details\n");
   printf("\n");
   printf("\n");
   printf("\n");
   fflush(stdout);

   FILE *fid2 = fopen("GTS_VERSION_INFO.out","wt");
   fprintf(fid2,"============ GIT INFORMATION ============\n");
   fprintf(fid2,"\n");
   fprintf(fid2,"       VERSION: %s\n", GIT_VERSION);
   fprintf(fid2,"        BRANCH: %s\n", GIT_BRANCH);
   fprintf(fid2,"    COMMIT MSG: %s\n", GIT_COM_MSG);
   fprintf(fid2,"    COMMIT SHA: %s\n", GIT_COM_SHA);
   fprintf(fid2,"      GIT DATE: %s\n", GIT_DATE);
   fprintf(fid2,"    BUILD DATE: %s\n", BUILD_DATE);
   fprintf(fid2,"  WORKING TREE: %s\n", GIT_WTREE);
   fprintf(fid2,"\n");
   fprintf(fid2,"=========================================\n");
   fprintf(fid2,"\n");
   fprintf(fid2,"If \"WORKING TREE\" is not clean, check \"GTS_GIT_DIFF.out\" for details\n");
   fprintf(fid2,"\n");
   fprintf(fid2,"\n");
   fprintf(fid2,"\n");
   fclose(fid2);
#endif

   FILE *fid= fopen("GTS_GIT_DIFF.out","w+");
   fprintf(fid,"%s\n",gts_git_diff_out);
   fclose(fid);

   return(0);
}
