#include "utility.h"

#undef __FUNCT__  
#define __FUNCT__ "nerror"
void nerror(char error_text[])
{
  fprintf(stderr, "Run time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  PetscFinalize();
  exit(1);
}

/* routine to allocate memory for 1d and 2d arrays */

#undef __FUNCT__  
#define __FUNCT__ "iVector1D"
int *iVector1D(int low, int high)
{
  /* allocate an integer vector with subscript [low,high] */ 
  int *p;
  if (high < low) high=low;
  p = (int *)malloc((unsigned) ((high-low+1)*sizeof(int)));
  if(!p) nerror("Allocation failure in iVector1D()");
  return p-low;
}

#undef __FUNCT__  
#define __FUNCT__ "iVector2D"
int **iVector2D(int row_low, int row_high, int column_low, int column_high)
{
  /*** 
     allocate an two-domensional array with subscript 
     [row_low,row_high][column_low, column_high] 
  ***/
  
  int i, nrow=row_high-row_low+1, ncolumn=column_high-column_low+1;
  int **p;

  /* allocate pointers to rows */
  p = (int **) malloc((unsigned)((nrow)*sizeof(int *)));
  if (!p) nerror("allocation failure 1 in iVector2D()");
  p -= row_low;

  /* allocate rows and set pointers to them */
  for(i=row_low;i<=row_high;i++) {
    p[i]=(int *)malloc((unsigned) (ncolumn)*sizeof(int));
    if (!p[i]) nerror("allocation failure 2 in iVector2D()");
    p[i] -= column_low;
  }
  return p;
}


#undef __FUNCT__
#define __FUNCT__ "iVector3D"
int ***iVector3D(int row_low, int row_high, int col_low, int col_high, int hei_low, int hei_high)
{
  /***
     allocate an three-domensional array with subscript
     [row_low,row_high][col_low, col_high][hei_low, hei_high]
  ***/
 
  int i, j;
  int nrow=row_high-row_low+1, ncol=col_high-col_low+1, nhei=hei_high-hei_low+1;
  int ***p;

  /* allocate pointers to rows */
  p = (int ***) malloc((unsigned)((nrow)*sizeof(int **)));
  if (!p) nerror("allocation failure 1 in iVector3D()");
  p -= row_low;

  /* allocate rows and set pointers to them */
  for(i=row_low;i<=row_high;i++) {
    p[i]=(int **)malloc((unsigned)((ncol)*sizeof(int *)));
    if (!p[i]) nerror("allocation failure 2 in iVector3D()");
    p[i] -= col_low;
  }

  /* allocate hei and set pointers to them */
  for(i=row_low;i<=row_high;i++) {           
     for(j=col_low;j<=col_high;j++) {           
        p[i][j]=(int *)malloc((unsigned)((nhei)*sizeof(int)));
        if (!p[i][j]) nerror("allocation failure 3 in iVector3D()");
        p[i][j] -= hei_low;            
     }
  }

  return p;
}


#undef __FUNCT__
#define __FUNCT__ "iVector4D"
int ****iVector4D(int row_low, int row_high, int col_low, int col_high, int hei_low, int hei_high, int hej_low, int hej_high)
{
  /***
     allocate an three-domensional array with subscript
     [row_low,row_high][col_low, col_high][hei_low, hei_high]
  ***/
 
  int i, j, k;
  int nrow=row_high-row_low+1, ncol=col_high-col_low+1, nhei=hei_high-hei_low+1, nhej=hej_high-hej_low+1;
  int ****p;

  /* allocate pointers to rows */
  p = (int ****) malloc((unsigned)((nrow)*sizeof(int ***)));
  if (!p) nerror("allocation failure 1 in iVector4D()");
  p -= row_low;

  /* allocate rows and set pointers to them */
  for(i=row_low;i<=row_high;i++) {
    p[i]=(int ***)malloc((unsigned) (ncol)*sizeof(int **));
    if (!p[i]) nerror("allocation failure 2 in iVector4D()");
    p[i] -= col_low;
  }

  /* allocate hei and set pointers to them */
  for(i=row_low;i<=row_high;i++) {           
     for(j=col_low;j<=col_high;j++) {           
        p[i][j]=(int **)malloc((unsigned) (nhei)*sizeof(int *));
        if (!p[i][j]) nerror("allocation failure 3 in iVector4D()");
        p[i][j] -= hei_low;            
     }
  }

  /* allocate hej and set pointers to them */
  for(i=row_low;i<=row_high;i++) {           
     for(j=col_low;j<=col_high;j++) {           
        for(k=hei_low;k<=hei_high;k++) {           
        p[i][j][k]=(int *)malloc((unsigned) (nhej)*sizeof(int ));
        if (!p[i][j][k]) nerror("allocation failure 4 in iVector4D()");
        p[i][j][k] -= hej_low;            
        }
     }
  }

  return p;
}



#undef __FUNCT__  
#define __FUNCT__ "dVector1D"
double *dVector1D(int nl, int nh)
{
  /* allocate a one dimensional double array with subscript [nl,nh] */
  double *v;
  
  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) nerror("allocation failure in dVector1D()");
  return v-nl;
}

#undef __FUNCT__  
#define __FUNCT__ "dVector2D"
double **dVector2D(int nrl,int nrh,int ncl,int nch)
{
  /* allocate a two dimensional double array with subscript [nrl,nrh][ncl,nch] */
  int i;
  double **m;
  
  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m) nerror("allocation failure 1 in dVector2D()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) nerror("allocation failure 2 in dVector2D()");
    m[i] -= ncl;
  }
  return m;
}

#undef __FUNCT__  
#define __FUNCT__ "dVector3D"
double ***dVector3D(int nrl, int nrh, int ncl, int nch, int nclrl, int nclrh)
{
  int     i, j;
  double  ***m;

  m=(double ***) malloc((unsigned) (nrh-nrl+1)*sizeof(double**));
  if (!m) nerror("allocation failure 1 in dVector3D()");
  m -= nrl;

  for(i=nrl; i<=nrh; i++) {
    m[i]=(double **) malloc((unsigned) (nch-ncl+1)*sizeof(double*));
    if (!m[i]) nerror("allocation failure 2 in dVector3D()");
    m[i] -= ncl;
  }

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      m[i][j]=(double *) malloc((unsigned) (nclrh-nclrl+1)*sizeof(double));
      if (!m[i][j]) nerror("allocation failure 3 in dVector3D()");
      m[i][j] -= nclrl;
    }
  }
  return m;
}

#undef __FUNCT__  
#define __FUNCT__ "dVector4D"
double ****dVector4D(int nrl, int nrh, int ncl, int nch, int nclrl, int nclrh, int nhl, int nhh)
{
  int     i, j, k;
  double  ****m;

  m=(double ****) malloc((unsigned) (nrh-nrl+1)*sizeof(double***));
  if (!m) nerror("allocation failure 1 in dVector4D()");
  m -= nrl;

  for(i=nrl; i<=nrh; i++) {
    m[i]=(double ***) malloc((unsigned) (nch-ncl+1)*sizeof(double**));
    if (!m[i]) nerror("allocation failure 2 in dVector4D()");
    m[i] -= ncl;
  }

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      m[i][j]=(double **) malloc((unsigned) (nclrh-nclrl+1)*sizeof(double*));
      if (!m[i][j]) nerror("allocation failure 3 in dVector4D()");
      m[i][j] -= nclrl;
    }
  }

  for(i=nrl; i<=nrh; i++) {
    for(j=ncl; j<=nch; j++) {
      for(k=nclrl; k<=nclrh; k++) {
        m[i][j][k]=(double *) malloc((unsigned) (nhh-nhl+1)*sizeof(double));
        if (!m[i][j][k]) nerror("allocation failure 4 in dVector4D()");
        m[i][j][k] -= nhl;
      }
    }
  }

  return m;
}

#undef __FUNCT__  
#define __FUNCT__ "MPI_Comm_vector"
MPI_Comm *MPI_Comm_vector(int low, int high)
{
  /* allocate an MPI_Comm array [low, high] */
  MPI_Comm *p;
  p = (MPI_Comm *)malloc((unsigned)((high-low+1)*sizeof(MPI_Comm)));
  if(!p) nerror("Allocation failure in iVector1D()");
  return p-low;  
}

#undef __FUNCT__  
#define __FUNCT__ "VecArray"
Vec *VecArray(int low, int high)
{
  /* allocate an Vec array [low, high] */
  Vec *p;
  p = (Vec *)malloc((unsigned)((high-low+1)*sizeof(Vec)));
  if(!p) nerror("Allocation failure in iVector1D()");
  return p-low;
}

#undef __FUNCT__  
#define __FUNCT__ "free_iVector1D"
void free_iVector1D(int *v, int nl, int nh)
{
	free((char*) (v+nl));
}

#undef __FUNCT__  
#define __FUNCT__ "free_iVector2D"
void free_iVector2D(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

#undef __FUNCT__  
#define __FUNCT__ "free_iVector3D"
void free_iVector3D(int ***m, int nrl, int nrh, int ncl, int nch, int nhil, int nhih)
{
	int i, j;

	for(i=nrh;i>=nrl;i--)
	  for(j=nch;i>=ncl;j--)
            free((char*) (m[i][j]+nhil));

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));

	free((char*) (m+nrl));
}

#undef __FUNCT__  
#define __FUNCT__ "free_iVector4D"
void free_iVector4D(int ****m, int nrl, int nrh, int ncl, int nch, int nhil, int nhih, int nhjl, int nhjh)
{
	int i, j, k;

	for(i=nrh;i>=nrl;i--)
	  for(j=nch;i>=ncl;j--)
	    for(k=nhih;k>=nhil;k--)
              free((char*) (m[i][j][k]+nhjl));

	for(i=nrh;i>=nrl;i--)
	  for(j=nch;i>=ncl;j--)
            free((char*) (m[i][j]+nhil));

	for(i=nrh;i>=nrl;i--) 
          free((char*) (m[i]+ncl));

	free((char*) (m+nrl));
}

#undef __FUNCT__  
#define __FUNCT__ "free_dVector1D"
void free_dVector1D(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}

#undef __FUNCT__  
#define __FUNCT__ "free_dVector2D"
void free_dVector2D(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

#undef __FUNCT__  
#define __FUNCT__ "free_dVector3D"
void free_dVector3D(double ***m, int nrl, int nrh, int ncl, int nch, int nclrl, int nclrh)
{
  int i, j;

  for (i=nrh; i>=nrl; i--)
    for (j=nch; j>=ncl; j--)
      free((char*) (m[i][j]+nclrl));

  for (i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl));

  free((char*) (m+nrl));
}

#undef __FUNCT__  
#define __FUNCT__ "free_dVector4D"
void free_dVector4D(double ****m, int nrl, int nrh, int ncl, int nch, int nclrl, int nclrh, int nhl, int nhh)
{
  int i, j, k;

  for (i=nrh; i>=nrl; i--)
    for (j=nch; j>=ncl; j--)
      for (k=nhh; j>=nhl; j--)
        free((char*) (m[i][j][k]+nhl));

  for (i=nrh; i>=nrl; i--)
    for (j=nch; j>=ncl; j--)
      free((char*) (m[i][j]+nclrl));

  for (i=nrh; i>=nrl; i--) 
    free((char*) (m[i]+ncl));

  free((char*) (m+nrl));
}

#undef __FUNCT__  
#define __FUNCT__ "free_Comm_vector"
void free_MPI_Comm_vector(MPI_Comm *p, int low, int high)
{
  int i;

  for (i=low; i<=high; i++)
    MPI_Comm_free(&p[i]);
}

#undef __FUNCT__  
#define __FUNCT__ "free_VecArray"
int free_VecArray(Vec *p, int low, int high)
{
  int i, ierr;

  for (i=low; i<=high; i++) {
    ierr = VecDestroy(&(p[i])); CHKERRQ(ierr);
  }
  return 0;
}
