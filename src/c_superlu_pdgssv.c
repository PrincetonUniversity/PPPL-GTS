/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include "pdsp_defs.h"
#include "util.h"
#include "Cnames.h"
#include <omp.h>
#include <math.h>

/* kind of integer to hold a pointer.  Use 'long int'
   so it works on 64-bit systems. */
typedef long int fptr;  /* 64 bit */

typedef struct {
    SuperMatrix *X;
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;

void
#ifdef _F_TO_C_NO_UNDERSCORE
c_superlu_pdgssv(int *mype,
#else
c_superlu_pdgssv_(int *mype,
#endif
                  int *nprocs, int *n, int *nnz, int *nrhs, double *values,
		 int *rowind, int *colptr, double *b, double *x, int *ldb, 
                 fptr *f_factors, /* a handle containing the address
                                     pointing to the factored matrices */
                 int *info,
                 int *iopt)
{
/*  static int iopt=0;*/

    SuperMatrix *A, *X, *L, *U;
    SCformat *Lstore;
    NCformat *Ustore;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      relax, panel_size, permc_spec, i, j;
    double   rpg, rcond;
    pdgstrf_options_t  pdgstrf_options;
    superlu_memusage_t superlu_memusage;
    factors_t *LUfactors;

    DNformat  *Xstore;
    double    *Xmat;

    trans_t   trant;
    Gstat_t   Gstat;
    double    t0;      /* temporary time */
    double    *utime;
    flops_t   *ops, flopcnt;



if(*iopt == 0) {
    /* Adjust to 0-based indexing */
#pragma omp parallel for private(i)
    for (i = 0; i < *nnz; ++i) --rowind[i];
#pragma omp parallel for private(i)
    for (i = 0; i <= *n; ++i) --colptr[i];

/*the matrix passed is not correct somehow, so check here -- start 
        if (*mype == 0) {
           FILE *fp;
           int i;
           char dumpname[100];
           sprintf(dumpname, "super_csc_0_%d", *mype);
           fp = fopen(dumpname, "w");
               fprintf(fp," nrow=%d \t ncol=%d \t nnz=%d \n", *n, *n, *nnz);
           for (i=0; i < *nnz; ++i)
               fprintf(fp," %d\t%d\t%e \n", i, rowind[i], values[i]);
           fclose(fp);
        }
 nrow=40107      ncol=40107      nnz=1017487 
 1017486        40106   2.044417e+00 
right on 10/23/06 --*/
/*the matrix passed is not correct somehow, so check here --  end  --*/

        A = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
        X = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
        L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
        U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    dCreate_CompCol_Matrix(A, *n, *n, *nnz, values, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(X, *n, *nrhs, x, *ldb, SLU_DN, SLU_D, SLU_GE);
fprintf(stderr, "\t%d_iopt=%d: _nrow=%d \t _ncol=%d \t _nnz=%d \t _f_factors=%ld\n",
        *mype, *iopt, A->nrow, A->ncol, *nnz, *f_factors);

    if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */    	
    permc_spec = 1;
    get_perm_c(permc_spec, A, perm_c);
fprintf(stderr, "\t%d_iopt=%d: get_perm_c\n", *mype, *iopt);

    pdgstrf_options.nprocs = *nprocs;
    pdgstrf_options.fact = EQUILIBRATE;
    pdgstrf_options.trans = NOTRANS;
    pdgstrf_options.refact = NO;
    pdgstrf_options.panel_size = sp_ienv(1);
    pdgstrf_options.relax = sp_ienv(2);
    pdgstrf_options.diag_pivot_thresh = 1.;
    pdgstrf_options.usepr = NO;
    pdgstrf_options.drop_tol = 0.;
    pdgstrf_options.perm_c = perm_c;
    pdgstrf_options.perm_r = perm_r;
    pdgstrf_options.work = NULL;
    pdgstrf_options.lwork = 0;
fprintf(stderr, "\t%d_iopt=%d: initial setup\n", *mype, *iopt);


/*
    pdgssv(*nprocs, &A, perm_c, perm_r, &L, &U, &B, info);
    c_pdgssvx(*mype, *nprocs, &pdgstrf_options, A, perm_c, perm_r, &equed, R, C, L, U, B, X, &rpg, &rcond, ferr, berr, &superlu_memusage, info);
*/
{/* start of solver */
    *info = 0;
    SuperMatrix AC; /* Matrix postmultiplied by Pc */

    /* ------------------------------------------------------------
       Allocate storage and initialize statistics variables. 
       ------------------------------------------------------------*/
    panel_size = pdgstrf_options.panel_size;
    relax = pdgstrf_options.relax;
    StatAlloc(*n, *nprocs, panel_size, relax, &Gstat);
    StatInit(*n, *nprocs, &Gstat);
    utime = Gstat.utime;
    ops = Gstat.ops;
fprintf(stderr, "\t%d_iopt=%d: initial solver\n", *mype, *iopt);
    
    /* ------------------------------------------------------------
       Convert A to NC format when necessary.
       ------------------------------------------------------------*/
    /* A->Stype == NC */
        trant = pdgstrf_options.trans;

    /* ------------------------------------------------------------
       Perform the LU factorization.
       ------------------------------------------------------------*/
    {   /* Obtain column etree, the column count (colcnt_h) and supernode
	   partition (part_super_h) for the Householder matrix. */
	t0 = SuperLU_timer_();
	sp_colorder(A, perm_c, &pdgstrf_options, &AC);
	utime[ETREE] = SuperLU_timer_() - t0;
        Destroy_SuperMatrix_Store(A);
fprintf(stderr, "\t%d_iopt=%d: A is destroyed to release memory\n", *mype, *iopt);

	/* Compute the LU factorization of A*Pc. */
	t0 = SuperLU_timer_();
	pdgstrf(&pdgstrf_options, &AC, perm_r, L, U, &Gstat, info);
	utime[FACT] = SuperLU_timer_() - t0;
fprintf(stderr, "\t%d_iopt=%d: L U formed\n", *mype, *iopt);
	
	flopcnt = 0;
	for (i = 0; i < *nprocs; ++i) flopcnt += Gstat.procstat[i].fcops;
	ops[FACT] = flopcnt;
fprintf(stderr, "	LU time=%8.2f \tflops=%e\tMflops=%8.2f\n", utime[FACT], ops[FACT], ops[FACT]*1e-6/utime[FACT]);
    }

    if ( *info > 0 ) {
        /* Compute the reciprocal pivot growth factor of the leading
	       rank-deficient *info columns of A. */
        exit(1);
    } else {
       Xstore = X->Store;
       Xmat   = Xstore->nzval;

	/* ------------------------------------------------------------
	   Compute the solution matrix X.
	   ------------------------------------------------------------*/
	for (j = 0; j < *nrhs; j++) {  /* Save a copy of the right hand sides */
#pragma omp parallel for private(i)
	    for (i = 0; i < *ldb; i++)
		Xmat[i + j*(*ldb)] = b[i + j*(*ldb)];
	}
    
	t0 = SuperLU_timer_();
	dgstrs(trant, L, U, perm_r, perm_c, X, &Gstat, info);
	utime[SOLVE] = SuperLU_timer_() - t0;
	ops[SOLVE] = ops[TRISOLVE];
        if (*mype == 0) fprintf(stderr,"	X time=%8.2f \tflops=%e\tMflops=%8.2f\n", utime[SOLVE], ops[SOLVE], ops[SOLVE]*1e-6/utime[SOLVE]);
    }
fprintf(stderr, "\t%d_iopt=%d: cleaned up\n", *mype, *iopt);

    /* ------------------------------------------------------------
       Deallocate storage after factorization.
       ------------------------------------------------------------*/
        SUPERLU_FREE(pdgstrf_options.etree);
        SUPERLU_FREE(pdgstrf_options.colcnt_h);
	SUPERLU_FREE(pdgstrf_options.part_super_h);
        Destroy_CompCol_Permuted(&AC);

    /* ------------------------------------------------------------
       Print timings, then deallocate statistic variables.
    PrintStat(&Gstat);
       ------------------------------------------------------------*/
    StatFree(&Gstat);
}/* end of solver*/

        /* Save the LU factors in the factors handle */
        LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
        LUfactors->X = X;
        LUfactors->L = L;
        LUfactors->U = U;
        LUfactors->perm_c = perm_c;
        LUfactors->perm_r = perm_r;
        *f_factors = (fptr) LUfactors;

    /* Restore to 1-based indexing */
#pragma omp parallel for private(i)
    for (i = 0; i < *nnz; ++i) ++rowind[i];
#pragma omp parallel for private(i)
    for (i = 0; i <= *n; ++i) ++colptr[i];



} else if(*iopt != 0) {
/* Extract the LU factors in the factors handle */
        LUfactors = (factors_t*) *f_factors;
        X = LUfactors->X;          /*1st time is correct*/
        L = LUfactors->L;          /*1st time is correct*/
        U = LUfactors->U;          /*1st time is correct*/
        perm_c = LUfactors->perm_c;
        perm_r = LUfactors->perm_r;
        if (*mype == 0)
fprintf(stderr,"_iopt=%d: _nrow=%d \t _ncol=%d \t _nnz=%d \t _f_factors=%ld\n", *iopt, L->nrow, L->ncol, *nnz, *f_factors);

    pdgstrf_options.nprocs = *nprocs;
    pdgstrf_options.fact = FACTORED;
    pdgstrf_options.trans = NOTRANS;
    pdgstrf_options.refact = YES;
    pdgstrf_options.panel_size = sp_ienv(1);
    pdgstrf_options.relax = sp_ienv(2);
    pdgstrf_options.diag_pivot_thresh = 1.;
    pdgstrf_options.usepr = YES;
    pdgstrf_options.drop_tol = 0.;
    pdgstrf_options.perm_c = perm_c;
    pdgstrf_options.perm_r = perm_r;
    pdgstrf_options.work = NULL;
    pdgstrf_options.lwork = 0;


/*  
    c_pdgssvx(*mype, *nprocs, pdgstrf_options_tmp, A, perm_c, perm_r, equed_tmp, R, C, L, U, B, X, &rpg, &rcond, ferr, berr, &superlu_memusage, info);
*/
{/*solver start*/
    *info = 0;

    /* ------------------------------------------------------------
       Allocate storage and initialize statistics variables. 
       ------------------------------------------------------------*/
    panel_size = pdgstrf_options.panel_size;
    relax = pdgstrf_options.relax;
    StatAlloc(*n, *nprocs, panel_size, relax, &Gstat);
    StatInit(*n, *nprocs, &Gstat);
    utime = Gstat.utime;
    ops = Gstat.ops;
    
    /* A->Stype == NC */
        trant = pdgstrf_options.trans;

    {
       Xstore = X->Store;
       Xmat   = Xstore->nzval;

	/* ------------------------------------------------------------
	   Compute the solution matrix X.
	   ------------------------------------------------------------*/
	for (j = 0; j < *nrhs; j++) {  /* Save a copy of the right hand sides */
#pragma omp parallel for private(i)
	    for (i = 0; i < *ldb; i++)
		Xmat[i + j*(*ldb)] = b[i + j*(*ldb)];
	}
    
	t0 = SuperLU_timer_();
	dgstrs(trant, L, U, perm_r, perm_c, X, &Gstat, info);
	utime[SOLVE] = SuperLU_timer_() - t0;
	ops[SOLVE] = ops[TRISOLVE];
        if (*mype == 0) fprintf(stderr,"	X time=%8.2f \tflops=%e\tMflops=%8.2f\n", utime[SOLVE], ops[SOLVE], ops[SOLVE]*1e-6/utime[SOLVE]);
    }

    /* ------------------------------------------------------------
       Print timings, then deallocate statistic variables.
    PrintStat(&Gstat);
       ------------------------------------------------------------*/
    StatFree(&Gstat);
}/*solver end*/

}/* end if(*iopt==1) {*/


    if ( *info == 0 ) {
#ifdef DEBUG_NOW
	Lstore = (SCformat *) L->Store;
	Ustore = (NCformat *) U->Store;
    	fprintf(stderr,"	#NZ in factor L = %d, #NZ in factor U = %d\n", Lstore->nnz, Ustore->nnz);
    	fprintf(stderr,"	L.row=%d, L.col=%d, U.row=%d, U.col = %d\n", L->nrow, L->ncol, U->nrow, U->ncol);
#endif
	
	superlu_QuerySpace(*nprocs, L, U, panel_size, &superlu_memusage);
	fprintf(stderr,"	L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	       superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6, superlu_memusage.expansions);
    } else {
        printf("dgssv() error returns INFO= %d\n", *info);
        exit(1);
    }

	for (j = 0; j < *nrhs; j++) {  /* Save a copy of the right hand sides */
#pragma omp parallel for private(i)
	    for (i = 0; i < *ldb; i++)
		x[i + j*(*ldb)] = Xmat[i + j*(*ldb)];
	}

#ifdef DEBUG_NOW
    double   rmsx, rmsb;
        c_rms(ldb, nrhs, x, &rmsx);
        c_rms(ldb, nrhs, b, &rmsb);
        if (*mype == 0) fprintf(stderr,"        iopt=%d, rmsb= %e, rmsx = %e\n", *iopt, rmsb, rmsx);
#endif

     fflush(stdout);
     fflush(stderr);
    /* *iopt++;
    return 0;*/
}


/* ---------------------------------------
In C, the first four elements of an array declared as X[3][3] are:
  X[0][0] X[0][1] X[0][2] X[1][0]

In Fortran, the first four elements are:
  X(1,1) X(2,1) X(3,1) X(1,2)

The order of indexing extends to any number of dimensions you declare. For example, the C declaration:
int arr1[2][10][15][20];

is equivalent to the Fortran declaration:
INTEGER arr1( 20, 15, 10, 2 )
--------------------------------------- */


int c_rms(int *ldb, int *nrhs, double *array, double *rms)
{
  int i, j;
  double total=0.;

        for (j = 0; j < *nrhs; j++) {  /* Save a copy of the right hand sides */
#pragma omp parallel for private(i) reduction(+:total)
            for (i = 0; i < *ldb; i++)
                 total += array[i + j*(*ldb)] * array[i + j*(*ldb)];
        }

  *rms = sqrt(total)/sqrt(*ldb);

  return 0;
}

