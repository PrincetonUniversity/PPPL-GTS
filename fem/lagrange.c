#include <petscksp.h>
#include "utility.h"


int lagrange_element_define_(int *order,
                             int *node,
                             int *ns,
                             int *nw
                               )
{
/*int node, order, nw, ns*/;
  int medge, nedge, nbubble;

/*  3 vertices defining a triangle 
    3 barycentric coordinates
*/
    *node = 3;

/*  number of collocation point for order 1 element
*/
    *order = 1;
    *nw =  3;
    *ns =  3;

/* higher element related
*/
    medge = *order - 1;   /* =1 for =2 =3 */
    nedge = medge * 3;  /* =3 for =6  =9 */
    nbubble = *ns - 3 - nedge;/* =0 for =1 =3 */
    printf("\tfinite element: order=%d ns=%d nw=%d medge=%d nedge=%d nbubble=%d \n",
           *order, *ns, *nw, medge,nedge,nbubble);
  return 0;
}


int lagrange_element_construct_(int *order,
                                int *ns,
                                int *nw,
                                double *weight,
                                double *NN,
                                double *dN1,
                                double *Eij
                               )
{
  int node /*, order, nw, ns*/;
  int medge, nedge, nbubble;
  int i, j, l, m, n;

/*  3 vertices defining a triangle 
    3 barycentric coordinates
    node = 3;
*/

/*  number of collocation point for order 1 element
    order = 1;
    *nw =  3;
    *ns =  3;
*/

/* higher element related
    medge = *order - 1;
    nedge = medge * 3;
    nbubble = *ns - 3 - nedge;
*/
    printf("\tfinitelement: lagrange_element_construct\n");

/* unit triangle related, never changes over each trianlge
    double phi[3][3];
    weight   = dVector1D(                   0, *nw -1);
    NN  = dVector2D(             0, *ns -1, 0, *nw -1);
    dN1 = dVector3D(       0, 2, 0, *ns -1, 0, *nw -1);
    Eij = dVector3D(  0, *ns -1, 0, *ns -1, 0, *nw -1);
*/

/*  phi(node, nw): phi_1 phi_2 phi_3 are valued barycentric coordinates 
                   at nw collocation points
*/
    double phi[3][3] = { 2./3. , 1./6. , 1./6.,
                         1./6. , 2./3. , 1./6.,
                         1./6. , 1./6. , 2./3.};
/*  weights at nw collocation points
*/
    weight[0] =1./3.;
    weight[1] =1./3.;
    weight[2] =1./3.;

/* Lagrange base functions evaluated at (ns, nw)
            NN[i][n] = phi[i][n];
*/
    for (n = 0; n < *nw; n++)
        for (i = 0; i < 3; i++)
            NN[i+n*3] = phi[i][n];

/* Lagrange base functions first derivative evaluated at (ns, nw)
    dN_i/dphi_l [n]
                dN1[l][i][n] = 0.;
                dN1[l][l][n] = 0.;
*/
    for (n = 0; n < *nw; n++)
        for (i = 0; i < *ns; i++)
            for (l = 0; l < 3; l++)
                dN1[l+ i*3 + n*(*ns)*3] = 0.;
    for (n = 0; n < *nw; n++)
        for (l = 0; l < 3; l++)
            dN1[l+l*3+n*3*(*ns)] = 1.;

/*  epsilon(3, 3) : epsilon(l,m) = [ lambda_l , lambda_m ]
*/
    double epsilon[3][3] = { 0., 1., -1., -1., 0., 1., 1., -1., 0. };

/* poisson bracket tensor
   ebrack___ = ___>___( d N_/ d lamda_)_( d N_/ d lamda_)_epsilon__
         ijn   l,m=1,3     i         l n     j         m n       lm
            Eij[i][j][n] = 0.;
*/
    for (n = 0; n < *nw; n++)
        for (j = 0; j < *ns; j++)
        for (i = 0; i < *ns; i++) {
            Eij[i + j*(*ns) + n*(*ns)*(*ns)] = 0.;
            for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
            Eij[i + j*(*ns) + n*(*ns)*(*ns)]  += 
            dN1[l + i*3 + n*(*ns)*3] * dN1[m + j*3 + n*(*ns)*3] * epsilon[l][m];
         }

  return 0;
}


int lagrange_element_matrix_(int *ns,
                             int *nw,
                             double *x,
                             double *y,
                             double *weight,
                             double *NN,
                             double *dN1,
                             double *Eij,
                             double *m,
                             double *s
                             )
{
/* matrix related, changes over each trianlge
*/
  int i, j, n;
/*double **m, **s;*/
  double *jacobian;
  double **dnx, **dny;

/*
  m      = dVector2D(0, *ns - 1, 0, *ns -1 );
  s      = dVector2D(0, *ns - 1, 0, *ns -1 );
*/
  jacobian   = dVector1D(          0, *nw-1);
  dnx        = dVector2D(0, *ns-1, 0, *nw-1);
  dny        = dVector2D(0, *ns-1, 0, *nw-1);

/* jacobian over each triangle
*/
  for (n = 0; n < *nw; ++n) {
      jacobian[n] = 0.;
      for (j = 0; j < *ns; ++j)
      for (i = 0; i < *ns; ++i)
          jacobian[n] += Eij[i + j*(*ns) + n*(*ns)*(*ns)] * x[i] * y[j];
  }

/* first order derivatives
 */
  for (n = 0; n < *nw; ++n)
      for (j = 0; j < *ns; ++j) {
          dnx[j][n] = 0.;
          dny[j][n] = 0.;
          for (i = 0; i < *ns; ++i) {
              dnx[j][n] += Eij[i + j*(*ns) + n*(*ns)*(*ns)] * y[i] / jacobian[n];
              dny[j][n] -= Eij[i + j*(*ns) + n*(*ns)*(*ns)] * x[i] / jacobian[n];
          }
      }

/* calculate mass (m), stiffness (s) matrices 
   over collocation points
          m[i][j]  = 0.;
          s[i][j]  = 0.;
*/
  for (j = 0; j < *ns; ++j)
  for (i = 0; i < *ns; ++i) {
      m[i + j*(*ns)]  = 0.;
      s[i + j*(*ns)]  = 0.;

      for (n = 0; n < *nw; ++n) {
          m[i + j*(*ns)] += weight[n] * jacobian[n] * NN[i + n*(*ns)] * NN[j + n*(*ns)];
          s[i + j*(*ns)] -= weight[n] * jacobian[n] * (dnx[i][n]*dnx[j][n] +  dny[i][n]*dny[j][n]);
      }
  }

  free_dVector1D(jacobian,          0, *nw-1);
  free_dVector2D(dnx, 0, *ns-1, 0, *nw-1);
  free_dVector2D(dny, 0, *ns-1, 0, *nw-1);
  return 0;
}
