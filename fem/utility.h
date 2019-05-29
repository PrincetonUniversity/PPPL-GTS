#ifndef UTILITY_H
#define UTILITY_H 
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <mpi.h>
#include "petscvec.h"

#define PI 3.141592653589793
 
/* memory allocation for 1d and 2d arrays */

void nerror(char error_text[]);

int *iVector1D(int, int);
int **iVector2D(int, int, int, int);
int ***iVector3D(int, int, int, int, int, int);
int ****iVector4D(int, int, int, int, int, int, int, int);
double *dVector1D(int, int);
double **dVector2D(int, int, int, int);
double ***dVector3D(int, int, int, int, int, int);
double ****dVector4D(int, int, int, int, int, int, int, int);
MPI_Comm *MPI_Comm_vector(int, int);
Vec *VecArray(int, int);

void free_iVector1D(int *, int, int);
void free_iVector2D(int **, int, int, int, int);
void free_iVector3D(int ***, int, int, int, int, int, int);
void free_iVector4D(int ****, int, int, int, int, int, int, int, int);
void free_dVector1D(double *, int, int);
void free_dVector2D(double **, int, int, int, int);
void free_dVector3D(double ***, int, int, int, int, int, int);
void free_dVector4D(double ****, int, int, int, int, int, int, int, int);
void free_MPI_Comm_vector(MPI_Comm *, int, int);
int free_VecArray(Vec *, int, int);

#endif

