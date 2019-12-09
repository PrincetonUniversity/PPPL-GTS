#include <stdio.h>
#include "gts_git_diff.h"
#include "gts_module_files.h"
#include <unistd.h>
#include <stdio.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <errno.h>
#include <netdb.h>
#include <time.h>

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



   // Print Out run-time information of GTS
   fid= fopen("GTS_RUNTIME_INFO.out","w+");


	// TIME information
	int hours, minutes, seconds, day, month, year;
	time_t now;
	time(&now);
	struct tm *local = localtime(&now);

   hours = local->tm_hour;      	// get hours since midnight (0-23)
   minutes = local->tm_min;     	// get minutes passed after the hour (0-59)
   seconds = local->tm_sec;     	// get seconds passed after minute (0-59)
   day = local->tm_mday;        	// get day of month (1 to 31)
   month = local->tm_mon + 1;   	// get month of year (0 to 11)
   year = local->tm_year + 1900;	// get year since 1900

	fprintf(fid,"Date: %02d/%02d/%d\n", day, month, year);
   fprintf(fid,"Time: %02d:%02d:%02d\n", hours, minutes, seconds);


	// Woring Directory information
   char cwd[PATH_MAX]; //current working directory
   if (getcwd(cwd, sizeof(cwd)) != NULL) {
      fprintf(fid,"Current working dir: %s\n", cwd);
   } else {
      perror("getcwd() error");
      return 1;
   }

   // hostname and IP information
   char hostbuffer[256];
   char *IPbuffer;
   struct hostent *host_entry;
   int hostname;
   // To retrieve hostname
   hostname = gethostname(hostbuffer, sizeof(hostbuffer));
   // To retrieve host information
   host_entry = gethostbyname(hostbuffer);
   // To convert an Internet network
   // address into ASCII string
   IPbuffer = inet_ntoa(*((struct in_addr*) host_entry->h_addr_list[0]));
   fprintf(fid,"Hostname: %s\n", hostbuffer);
   fprintf(fid,"Host IP: %s\n", IPbuffer);


   fclose(fid);

   return(0);
}
