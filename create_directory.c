#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

/* int mkdir(const char *pathname, mode_t mode); */

int create_directory(const char *pathname)
{
   int ierr;

   ierr = mkdir( pathname, 0700);
   if (ierr != 0)fprintf(stderr,"Error creating directory %s : ierr=%d  errno=%d   %s\n",pathname,ierr,errno,strerror(errno));

 /* if directory already exists, send back -2 to Fortran function call */
   if (errno == EEXIST)ierr=-2;

   return ierr;
}

int create_directory_(const char *pathname)
{
   int ierr;

   ierr = mkdir( pathname, 0700);
   if (ierr != 0)fprintf(stderr,"Error creating directory %s : ierr=%d  errno=%d   %s\n",pathname,ierr,errno,strerror(errno));

 /* if directory already exists, send back -2 to Fortran function call */
   if (errno == EEXIST)ierr=-2;

   return ierr;
}

