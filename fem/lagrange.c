#include <petscksp.h>
#include "utility.h"


/*int tric1_element_define_(int *order,
                          int *node,
                          int *nw,
                          int *ns,
                          int *ndof,
                          int *nvar*/
int lagrange_element_define_(int *order,
                          int *node,
                          int *ns,
                          int *nw,
                          int *ndof,
                          int *nvar
                         )
{
/*int node, order, nw, ns*/;
  int medge, nedge, nbubble;

/*3 vertices defines a triangle*/
    *node = 3;

/*number of collocation point for order 1 element */
    *nw =  25;

/*number of geometrical nodes*/
    if(*order == 2) {
       *ns =  6;
    } else {
       *ns =  3;
    }

/*number of dof on attached to each geometrical node*/
    if(*order > 2) {
       *ndof =  6;
    } else {
       *ndof =  1;
    }

/*number of dof on attached to each triangle*/
    *nvar =  (*ndof) * (*ns);

/*higher element related*/
    medge = 0; //*order - 1;   /* =1 for =2 =3 */
    nedge = 0; //medge * 3;  /* =3 for =6  =9 */
    nbubble = 0; //*ns - 3 - nedge;/* =0 for =1 =3 */
    PetscPrintf(PETSC_COMM_WORLD,
        "\tfinite element: order=%d node=%d nw=%d ns=%d ndof=%d nvar=%d\n",
        *order, *node, *nw, *ns, *ndof, *nvar);
  return 0;
}


/*int tric1_element_construct_(int *order,
                             int *node,
                             int *nw,
                             int *ns,
                             int *ndof,
                             int *nvar,
                             double *weight,
                             double *bNN,
                             double *dbN1,
                             double *dbN2,
                             double *bEij,*/
int lagrange_element_construct_(int *order,
                             int *ns,
                             int *nw,
                             int *ndof,
                             int *nvar,
                             double *weight,
                             double *NN,
                             double *dN1,
                             double *Eij,
                             double *bEij
                            )
{
  /*int node , order, nw, ns*/;
  int medge, nedge, nbubble;
  int i, j, l, m, n, li;
  PetscScalar ph0, ph1, ph2;

/*  3 vertices defining a triangle 
    3 barycentric coordinates
    node = 3;
*/

/*  number of collocation point for order 1 element
    order = 1;
    *nw =  25;
    *ns =  3;
*/

/* higher element related
    medge = *order - 1;
    nedge = medge * 3;
    nbubble = *ns - 3 - nedge;
*/
    PetscPrintf(PETSC_COMM_WORLD,"\tfinit element: tric1_element_construct\n");

/* unit triangle related, never changes over each trianlge
    double phi[3][3];
    weight   = dVector1D(                   0, *nw -1);
    NN  = dVector2D(             0, *ns -1, 0, *nw -1);
    dN1 = dVector3D(       0, 2, 0, *ns -1, 0, *nw -1);
    Eij = dVector3D(  0, *ns -1, 0, *ns -1, 0, *nw -1);
*/

/*  phi(node, nw): phi_1 phi_2 phi_3 are valued barycentric coordinates 
                   at nw collocation points
    double phi[3][3] = { 2./3. , 1./6. , 1./6.,
                         1./6. , 2./3. , 1./6.,
                         1./6. , 1./6. , 2./3.};
*/
/* 5th order - 25 points */
     PetscScalar a0,a1,a2,a3,a4,a5;
     PetscScalar b0,b1,b2,b3,b4,b5;
     PetscScalar c0,c1,c2,c3,c4,c5;
     a0 = 1./3.;
     a1 = 0.028844733232685; b1 = 0.485577633383657; c1 = 0.485577633383657;
     a2 = 0.781036849029926; b2 = 0.109481575485037; c2 = 0.109481575485037;
     a3 = 0.141707219414880; b3 = 0.307939838764121; c3 = 0.550352941820999;
     a4 = 0.025003534762686; b4 = 0.246672560639903; c4 = 0.728323904597411;
     a5 = 0.009540815400299; b5 = 0.066803251012200; c5 = 0.923655933587500;
PetscScalar phi[3][25]={a0,a1,b1,c1,a2,b2,c2,a3,b3,b3,c3,c3,a3,a4,b4,b4,c4,c4,a4,a5,b5,b5,c5,c5,a5,
                        a0,b1,c1,a1,b2,c2,a2,b3,a3,c3,b3,a3,c3,b4,a4,c4,b4,a4,c4,b5,a5,c5,b5,a5,c5,
                        a0,c1,a1,b1,c2,a2,b2,c3,c3,a3,a3,b3,b3,c4,c4,a4,a4,b4,b4,c5,c5,a5,a5,b5,b5};

/*  weights at nw collocation points
    weight[0] =1./3.;
    weight[1] =1./3.;
    weight[2] =1./3.;
*/
     PetscScalar w0,w1,w2,w3,w4,w5;
     w0 = 0.090817990382754;
     w1 = 0.036725957756467;
     w2 = 0.045321059435528;
     w3 = 0.072757916845420;
     w4 = 0.028327242531057;
     w5 = 0.009421666963733;
          n = 0   ;                   weight[n] = w0;
     for (n = 0+ 1; n < 0+ 1+3; n++)  weight[n] = w1;
     for (n = 0+ 4; n < 0+ 4+3; n++)  weight[n] = w2;
     for (n = 0+ 7; n < 0+ 7+6; n++)  weight[n] = w3;
     for (n = 0+13; n < 0+13+6; n++)  weight[n] = w4;
     for (n = 0+19; n < 0+19+6; n++)  weight[n] = w5;

/*epsilon(3, 3) : epsilon(l,m) = [ lambda_l , lambda_m ]*/
    double epsilon[3][3] = { 0., 1., -1., -1., 0., 1., 1., -1., 0. };

  PetscScalar *bNN;
  PetscScalar *dbN1;
    bNN = dVector1D(0, (*ns)*(*nw) - 1);
    dbN1 = dVector1D(0, (*ns)*(*ns)*(*nw) - 1);
/*Lagrange base functions as subparameter mapping evaluated at (ns, nw)
*/
if(*order == 2) {
/*
  <N>=< N1 N2 N3 N4 N5 N6 >
*/
    for (n = 0; n < *nw; n++) {
        ph0 = phi[0][n];
        ph1 = phi[1][n];
        ph2 = phi[2][n];

        bNN[0*(*nw)+n] = - ph0*(1. - 2.*ph0);
        bNN[1*(*nw)+n] =  4.*ph1*ph0;
        bNN[2*(*nw)+n] = - ph1*(1. - 2.*ph1);
        bNN[3*(*nw)+n] =  4.*ph1*ph2;
        bNN[4*(*nw)+n] = - ph2*(1. - 2.*ph2);
        bNN[5*(*nw)+n] =  4.*ph2*ph0;
       
        //0:dxi 1:deta
        dbN1[0*(*ns)*(*nw) + 0*(*nw) + n] = 1. - 4.*ph0;
        dbN1[1*(*ns)*(*nw) + 0*(*nw) + n] = 1. - 4.*ph0;

        dbN1[0*(*ns)*(*nw) + 1*(*nw) + n] = 4.*(ph0 - ph1);
        dbN1[1*(*ns)*(*nw) + 1*(*nw) + n] = - 4.*ph1;

        dbN1[0*(*ns)*(*nw) + 2*(*nw) + n] = -1. + 4.*ph1;
        dbN1[1*(*ns)*(*nw) + 2*(*nw) + n] = 0.;

        dbN1[0*(*ns)*(*nw) + 3*(*nw) + n] = 4.*ph2;
        dbN1[1*(*ns)*(*nw) + 3*(*nw) + n] = 4.*ph1;

        dbN1[0*(*ns)*(*nw) + 4*(*nw) + n] = 0.;
        dbN1[1*(*ns)*(*nw) + 4*(*nw) + n] = -1. + 4.*ph2;

        dbN1[0*(*ns)*(*nw) + 5*(*nw) + n] = -4.*ph2;
        dbN1[1*(*ns)*(*nw) + 5*(*nw) + n] = 4.*(ph0 - ph2);
       
        /*0:dxi^2 1:dxideta 2:deta^2
        dbN2[0*(*ns)*(*nw) + 0*(*nw) + n] = 4.;
        dbN2[1*(*ns)*(*nw) + 0*(*nw) + n] = 4.;
        dbN2[2*(*ns)*(*nw) + 0*(*nw) + n] = 4.;

        dbN2[0*(*ns)*(*nw) + 1*(*nw) + n] = -8.;
        dbN2[1*(*ns)*(*nw) + 1*(*nw) + n] = -4.;
        dbN2[2*(*ns)*(*nw) + 1*(*nw) + n] =  0.;

        dbN2[0*(*ns)*(*nw) + 2*(*nw) + n] = 4.;
        dbN2[1*(*ns)*(*nw) + 2*(*nw) + n] = 0.;
        dbN2[2*(*ns)*(*nw) + 2*(*nw) + n] = 0.;

        dbN2[0*(*ns)*(*nw) + 3*(*nw) + n] = 0.;
        dbN2[1*(*ns)*(*nw) + 3*(*nw) + n] = 4.;
        dbN2[2*(*ns)*(*nw) + 3*(*nw) + n] = 0.;

        dbN2[0*(*ns)*(*nw) + 4*(*nw) + n] = 0.;
        dbN2[1*(*ns)*(*nw) + 4*(*nw) + n] = 0.;
        dbN2[2*(*ns)*(*nw) + 4*(*nw) + n] = 4.;

        dbN2[0*(*ns)*(*nw) + 5*(*nw) + n] =  0.;
        dbN2[1*(*ns)*(*nw) + 5*(*nw) + n] = -4.;
        dbN2[2*(*ns)*(*nw) + 5*(*nw) + n] = -8.;*/
    }
} else {
/*
  <N>=<  1-ksi-eta   ksi   eta  >
*/
    for (i = 0; i < (*ns); i++)
        for (n = 0; n < *nw; n++)
            bNN[i*(*nw) + n] = phi[i][n];

    for (l = 0; l < 2; l++)
        for (i = 0; i < *ns; i++)
            for (n = 0; n < *nw; n++)
                dbN1[l*(*ns)*(*nw)+i*(*nw)+n] = 0.;

        i=0; // N1=1-ksi-eta
    for (n = 0; n < *nw; n++)
        dbN1[0*(*ns)*(*nw) + i*(*nw) + n] = -1.;
    for (n = 0; n < *nw; n++)
        dbN1[1*(*ns)*(*nw) + i*(*nw) + n] = -1.;

        i=1; // N2=ksi
    for (n = 0; n < *nw; n++)
        dbN1[0*(*ns)*(*nw) + i*(*nw) + n] =  1.;
    for (n = 0; n < *nw; n++)
        dbN1[1*(*ns)*(*nw) + i*(*nw) + n] =  0.;

        i=2; // N3=eta
    for (n = 0; n < *nw; n++)
        dbN1[0*(*ns)*(*nw) + i*(*nw) + n] =  0.;
    for (n = 0; n < *nw; n++)
        dbN1[1*(*ns)*(*nw) + i*(*nw) + n] =  1.;
}

/*poisson bracket tensor: Eij=dNi/dksi dNj/deta - dNi/deta dNj/dksi
   0   1  -1
  -1   0   1
   1  -1   0
*/
    for (i = 0; i < (*ns); i++)
    for (j = 0; j < (*ns); j++)
        for (n = 0; n < *nw; n++) {
            bEij[i*(*ns)*(*nw)+j*(*nw)+n] = dbN1[0*(*ns)*(*nw) + i*(*nw) + n] * dbN1[1*(*ns)*(*nw) + j*(*nw) + n]
                                          - dbN1[1*(*ns)*(*nw) + i*(*nw) + n] * dbN1[0*(*ns)*(*nw) + j*(*nw) + n];
        }

if(*order > 2) {
    /*gtsc1 base functions evaluated at (nvar, nw)*/
    for (n = 0; n < *nw; n++) {
        ph0 = phi[0][n];
        ph1 = phi[1][n];
        ph2 = phi[2][n];

        NN[ 0*(*nw)+n] = ph0*ph0*(10.*ph0 - 15.*ph0*ph0 + 6.*ph0*ph0*ph0 + 30.*ph1*ph2*(ph1+ph2));
        NN[ 1*(*nw)+n] = ph1*ph0*ph0*(3. - 2.*ph0 - 3.*ph1*ph1 + 6.*ph1*ph2);
        NN[ 2*(*nw)+n] = ph2*ph0*ph0*(3. - 2.*ph0 - 3.*ph2*ph2 + 6.*ph1*ph2);
        NN[ 3*(*nw)+n] = .5 * ph1*ph1*ph0*ph0*(1. - ph1 + 2.*ph2);
        NN[ 4*(*nw)+n] = ph1*ph2*ph0*ph0;
        NN[ 5*(*nw)+n] = .5 * ph2*ph2*ph0*ph0*(1. - ph2 + 2.*ph1);

        NN[ 6*(*nw)+n] = ph1*ph1*(10.*ph1 - 15.*ph1*ph1 + 6.*ph1*ph1*ph1 + 15.*ph2*ph2*ph0);
        NN[ 7*(*nw)+n] = .5 * ph1*ph1*(-8.*ph1 + 14.*ph1*ph1 - 6.*ph1*ph1*ph1 - 15.*ph2*ph2*ph0);
        NN[ 8*(*nw)+n] = .5 * ph1*ph1*ph2*(6. - 4.*ph1 - 3.*ph2 - 3.*ph2*ph2 + 3.*ph1*ph2);
        NN[ 9*(*nw)+n] = .25 * ph1*ph1*(2.*ph1*(1.-ph1)*(1.-ph1) + 5.*ph2*ph2*ph0);
        NN[10*(*nw)+n] = .5 * ph1*ph1*ph2*(-2. + 2.*ph1 + 1.*ph2 + ph2*ph2 - ph1*ph2);
        NN[11*(*nw)+n] = .25 * ph1*ph1*ph2*ph2*ph0 + .5 * ph1*ph1*ph1*ph2*ph2;

        NN[12*(*nw)+n] = ph2*ph2*(10.*ph2 - 15.*ph2*ph2 + 6.*ph2*ph2*ph2 + 15.*ph1*ph1*ph0);
        NN[13*(*nw)+n] = .5 * ph1*ph2*ph2*(6. - 3.*ph1 - 4.*ph2 - 3.*ph1*ph1 + 3.*ph1*ph2);
        NN[14*(*nw)+n] = .5 * ph2*ph2*(-8.*ph2 + 14.*ph2*ph2 - 6.*ph2*ph2*ph2 - 15.*ph1*ph1*ph0);
        NN[15*(*nw)+n] = .25 * ph1*ph1*ph2*ph2*ph0 + .5 * ph1*ph1*ph2*ph2*ph2;
        NN[16*(*nw)+n] = .5 * ph1*ph2*ph2*(-2. + 1.*ph1 + 2.*ph2 + ph1*ph1 - ph1*ph2);
        NN[17*(*nw)+n] = .25 * ph2*ph2*(2.*ph2*(1.-ph2)*(1.-ph2) + 5.*ph1*ph1*ph0);
    }

/*dNdphi at the collocation points phi(n, i), n=1:nw, i=1:nvar
  dN1( 0, i, n ) = d N_i / d ksi
  dN1( 1, i, n ) = d N_i / d eta
*/
    for (l = 0; l < 2; l++)
        for (i = 0; i < *nvar; i++)
            for (n = 0; n < *nw; n++)
                dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = 0.;

    for (n = 0; n < *nw; n++) {
        ph0 = phi[0][n];
        ph1 = phi[1][n];
        ph2 = phi[2][n];

        i=1-1; //N1
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = (-2.)*ph0*(10.*ph0 - 15.*ph0*ph0 + 6.*ph0*ph0*ph0 + 30.*ph1*ph2*(ph1+ph2))
                     + ph0*ph0*(10.*(-1.) - 15.*(-2.)*ph0 + 6.*(-3.)*ph0*ph0 + 30.*ph2*(ph1+ph2) + 30.*ph1*ph2);
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = (-2.)*ph0*(10.*ph0 - 15.*ph0*ph0 + 6.*ph0*ph0*ph0 + 30.*ph1*ph2*(ph1+ph2))
                     + ph0*ph0*(10.*(-1.) - 15.*(-2.)*ph0 + 6.*(-3.)*ph0*ph0 + 30.*ph1*(ph1+ph2) + 30.*ph1*ph2);

        i=2-1; //N2
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = ph0*ph0*(3. - 2.*ph0 - 3.*ph1*ph1 + 6.*ph1*ph2)
                     + ph1*(-2.)*ph0*(3. - 2.*ph0 - 3.*ph1*ph1 + 6.*ph1*ph2)
                     + ph1*ph0*ph0*(3*(0.) - 2.*(-1.) -3.*(2.)*ph1 + 6.*ph2);
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = ph1*(-2.)*ph0*(3. - 2.*ph0 - 3.*ph1*ph1 + 6.*ph1*ph2)
                     + ph1*ph0*ph0*(3*(0.) - 2.*(-1.) -3.*(0.) + 6.*ph1);

        i=3-1; //N3
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = ph2*(-2.)*ph0*(3. - 2.*ph0 - 3.*ph2*ph2 + 6.*ph1*ph2)
                     + ph2*ph0*ph0*(3.*(0.) - 2.*(-1.) -3.*(0.) + 6.*ph2);
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] =     ph0*ph0*(3. - 2.*ph0 - 3.*ph2*ph2 + 6.*ph1*ph2)
                     + ph2*(-2.)*ph0*(3. - 2.*ph0 - 3.*ph2*ph2 + 6.*ph1*ph2)
                     + ph2*ph0*ph0*(3.*(0.) - 2.*(-1.) - 3.*(2.)*ph2 + 6.*ph1);

        i=4-1; //N4
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*(2.)*ph1*ph0*ph0*(1. - ph1 + 2.*ph2)
                     + .5*ph1*ph1*(-2.)*ph0*(1. - ph1 + 2.*ph2)
                     + .5*ph1*ph1*ph0*ph0*(1.*(0.) - 1. + 2.*(0.));
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*ph1*ph1*(-2.)*ph0*(1. - ph1 + 2.*ph2)
                     + .5*ph1*ph1*ph0*ph0*(1.*(0.) - 1.*(0.) + 2.*(1.));

        i=5-1; //N5
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = ph2*ph0*ph0 + ph1*ph2*(-2.)*ph0;
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = ph1*ph0*ph0 + ph1*ph2*(-2.)*ph0;

        i=6-1; //N6
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*ph2*ph2*(-2.)*ph0*(1. + 2.*ph1 - ph2)
                     + .5*ph2*ph2*ph0*ph0*(1.*(0.) + 2.*(1.) - 0.);
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*(2.)*ph2*ph0*ph0*(1. + 2.*ph1 - ph2)
                     + .5*ph2*ph2*(-2.)*ph0*(1. + 2.*ph1 - ph2)
                     + .5*ph2*ph2*ph0*ph0*(1.*(0.) + 2.*(0.) - 1.);

        i=7-1; //N7
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = (2.)*ph1*(10.*ph1 - 15.*ph1*ph1 + 6.*ph1*ph1*ph1 + 15.*ph2*ph2*ph0)
                     + ph1*ph1*(10.*(1.) - 15.*(2.)*ph1 + 6.*(3.)*ph1*ph1 + 15.*ph2*ph2*(-1.));
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = ph1*ph1*(10.*(0.) - 15.*(0.) + 6.*(0.) + 15.*(2.)*ph2*ph0 + 15.*ph2*ph2*(-1.));

        i=8-1; //N8
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*(2.)*ph1*(-8.*ph1 + 14.*ph1*ph1 - 6.*ph1*ph1*ph1 - 15.*ph2*ph2*ph0)
                     + .5*ph1*ph1*(-8.*(1.) + 14.*(2.)*ph1 - 6.*(3.)*ph1*ph1 - 15.*ph2*ph2*(-1.));
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*ph1*ph1*(-8.*(0.) + 14.*(0.) - 6.*(0.) - 15.*(2.)*ph2*ph0 - 15.*ph2*ph2*(-1.));

        i=9-1; //N9
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*(2.)*ph1*ph2*(6. - 4.*ph1 - 3.*ph2 - 3.*ph2*ph2 + 3.*ph1*ph2)
                     + .5*ph1*ph1*ph2*(6.*(0.) - 4.*(1.) - 3.*(0.) - 3.*(0.) + 3.*ph2);
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*ph1*ph1*(1.)*(6. - 4.*ph1 - 3.*ph2 - 3.*ph2*ph2 + 3.*ph1*ph2)
                     + .5*ph1*ph1*ph2*(6.*(0.) - 4.*(0.) - 3.*1. - 3.*(2.)*ph2 + 3.*ph1);

        i=10-1; //N10
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .25*(2.)*ph1*(2.*ph1*(1.-ph1)*(1.-ph1) + 5.*ph2*ph2*ph0)
                     + .25*ph1*ph1*(2.*(1-ph1)*(1-ph1) + 2.*ph1*(-2.)*(1.-ph1) + 5.*ph2*ph2*(-1.));
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .25*ph1*ph1*(2.*(0.) + 5.*(2.)*ph2*ph0 + 5.*ph2*ph2*(-1.));

        i=11-1; //N11
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*(2.)*ph1*ph2*(-2. + 2.*ph1 + ph2 + ph2*ph2 - ph1*ph2)
                     + .5*ph1*ph1*ph2*(-2.*(0.) + 2.*(1.) + 0. + 0. - ph2);
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*ph1*ph1*(1.)*(-2. + 2.*ph1 + ph2 + ph2*ph2 - ph1*ph2)
                     + .5*ph1*ph1*ph2*(-2.*(0.) + 2.*(0.) + 1. + 2.*ph2 - ph1);

        i=12-1; //N12
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .25*(2.)*ph1*ph2*ph2*ph0 + .25*ph1*ph1*ph2*ph2*(-1.) + .5*(3.)*ph1*ph1*ph2*ph2;
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .25*ph1*ph1*(2.)*ph2*ph0 + .25*ph1*ph1*ph2*ph2*(-1.) + .5*ph1*ph1*ph1*(2.)*ph2;

        i=13-1; //N13
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = ph2*ph2 * (10.*(0.) - 15.*(0.) + 6.*(0.) + 15.*(2.)*ph1*ph0 + 15.*ph1*ph1*(-1.));
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = (2.)*ph2*(10.*ph2 - 15.*ph2*ph2 + 6.*ph2*ph2*ph2 + 15.*ph1*ph1*ph0)
                     + ph2*ph2*(10.*(1.) - 15.*(2.)*ph2 + 6.*(3.)*ph2*ph2 + 15.*ph1*ph1*(-1.));

        i=14-1; //N14
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*(1.)*ph2*ph2*(6. - 3.*ph1 - 4.*ph2 - 3.*ph1*ph1 + 3.*ph1*ph2)
                     + .5*ph1*ph2*ph2*(6.*(0.) - 3.*(1.) - 4.*(0.) - 3.*(2.)*ph1 + 3.*(1.)*ph2);
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*ph1*(2.)*ph2*(6. - 3.*ph1 - 4.*ph2 - 3.*ph1*ph1 + 3.*ph1*ph2)
                     + .5*ph1*ph2*ph2*(6.*(0.) - 3.*(0.) - 4.*(1.) - 3.*(0.) + 3.*ph1*(1.));

        i=15-1; //N15
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*ph2*ph2*(-8.*(0.) + 14.*(0.) - 6.*0. - 15.*(2.)*ph1*ph0 - 15.*ph1*ph1*(-1.));
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*(2.)*ph2*(-8.*ph2 + 14.*ph2*ph2 - 6.*ph2*ph2*ph2 - 15.*ph1*ph1*ph0)
                     + .5*ph2*ph2*(-8.*(1.) + 14.*(2.)*ph2 - 6.*(3.)*ph2*ph2 - 15.*ph1*ph1*(-1.));

        i=16-1; //N16
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .25*(2.)*ph1*ph2*ph2*ph0 + .25*ph1*ph1*ph2*ph2*(-1.) + .5*(2.)*ph1*ph2*ph2*ph2;
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .25*ph1*ph1*(2.)*ph2*ph0 + .25*ph1*ph1*ph2*ph2*(-1.) + .5*ph1*ph1*(3.)*ph2*ph2;

        i=17-1; //N17
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*(1.)*ph2*ph2*(-2. + ph1 + 2.*ph2 + ph1*ph1 - ph1*ph2) 
                     + .5*ph1*ph2*ph2*(-2.*(0.) + 1. + 2.*(0.) + (2.)*ph1 - ph2);
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .5*ph1*(2.)*ph2*(-2. + ph1 + 2.*ph2 + ph1*ph1 - ph1*ph2) 
                     + .5*ph1*ph2*ph2*(-2.*(0.) + 0. + 2.*(1.) + 0. - ph1);

        i=18-1; //N18
        l=0; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .25*ph2*ph2*(2.*0. + 5.*(2.)*ph1*ph0 + 5.*ph1*ph1*(-1.));
        l=1; li=l+i*2+n*2*(*nvar);
        dN1[l*(*nvar)*(*nw)+i*(*nw)+n] = .25*(2.)*ph2*(2.*ph2*(1.-ph2)*(1.-ph2) + 5.*ph1*ph1*ph0)
                     + .25*ph2*ph2*(2.*(1.-ph2)*(1.-ph2) + 2.*ph2*(-2.)*(1.-ph2) + 5.*ph1*ph1*(-1.));
    }
} else if(*order == 2) {
    for (n = 0; n < *nw; n++) {
        ph0 = phi[0][n];
        ph1 = phi[1][n];
        ph2 = phi[2][n];

        NN[0*(*nw)+n] = - ph0*(1. - 2.*ph0);
        NN[1*(*nw)+n] =  4.*ph1*ph0;
        NN[2*(*nw)+n] = - ph1*(1. - 2.*ph1);
        NN[3*(*nw)+n] =  4.*ph1*ph2;
        NN[4*(*nw)+n] = - ph2*(1. - 2.*ph2);
        NN[5*(*nw)+n] =  4.*ph2*ph0;

        //0:dxi 1:deta
        dN1[0*(*ns)*(*nw) + 0*(*nw) + n] = 1. - 4.*ph0;
        dN1[1*(*ns)*(*nw) + 0*(*nw) + n] = 1. - 4.*ph0;

        dN1[0*(*ns)*(*nw) + 1*(*nw) + n] = 4.*(ph0 - ph1);
        dN1[1*(*ns)*(*nw) + 1*(*nw) + n] = - 4.*ph1;

        dN1[0*(*ns)*(*nw) + 2*(*nw) + n] = -1. + 4.*ph1;
        dN1[1*(*ns)*(*nw) + 2*(*nw) + n] = 0.;

        dN1[0*(*ns)*(*nw) + 3*(*nw) + n] = 4.*ph2;
        dN1[1*(*ns)*(*nw) + 3*(*nw) + n] = 4.*ph1;

        dN1[0*(*ns)*(*nw) + 4*(*nw) + n] = 0.;
        dN1[1*(*ns)*(*nw) + 4*(*nw) + n] = -1. + 4.*ph2;

        dN1[0*(*ns)*(*nw) + 5*(*nw) + n] = -4.*ph2;
        dN1[1*(*ns)*(*nw) + 5*(*nw) + n] = 4.*(ph0 - ph2);

        /*0:dxi^2 1:dxideta 2:deta^2
        dN2[0*(*ns)*(*nw) + 0*(*nw) + n] = 4.;
        dN2[1*(*ns)*(*nw) + 0*(*nw) + n] = 4.;
        dN2[2*(*ns)*(*nw) + 0*(*nw) + n] = 4.;

        dN2[0*(*ns)*(*nw) + 1*(*nw) + n] = -8.;
        dN2[1*(*ns)*(*nw) + 1*(*nw) + n] = -4.;
        dN2[2*(*ns)*(*nw) + 1*(*nw) + n] =  0.;

        dN2[0*(*ns)*(*nw) + 2*(*nw) + n] = 4.;
        dN2[1*(*ns)*(*nw) + 2*(*nw) + n] = 0.;
        dN2[2*(*ns)*(*nw) + 2*(*nw) + n] = 0.;

        dN2[0*(*ns)*(*nw) + 3*(*nw) + n] = 0.;
        dN2[1*(*ns)*(*nw) + 3*(*nw) + n] = 4.;
        dN2[2*(*ns)*(*nw) + 3*(*nw) + n] = 0.;

        dN2[0*(*ns)*(*nw) + 4*(*nw) + n] = 0.;
        dN2[1*(*ns)*(*nw) + 4*(*nw) + n] = 0.;
        dN2[2*(*ns)*(*nw) + 4*(*nw) + n] = 4.;

        dN2[0*(*ns)*(*nw) + 5*(*nw) + n] =  0.;
        dN2[1*(*ns)*(*nw) + 5*(*nw) + n] = -4.;
        dN2[2*(*ns)*(*nw) + 5*(*nw) + n] = -8.;*/
    }
} else {
/*c1 base functions evaluated at (nvar, nw)*/
    for (n = 0; n < *nw; n++) {
        NN[ 0*(*nw)+n] = phi[0][n];
        NN[ 1*(*nw)+n] = phi[1][n];
        NN[ 2*(*nw)+n] = phi[2][n];
    }

/*dNdphi at the collocation points phi(n, i), n=1:nw, i=1:nvar
  dN1( 0, i, n ) = d N_i / d ksi
  dN1( 1, i, n ) = d N_i / d eta
*/
    for (l = 0; l < 2; l++)
        for (i = 0; i < *ns; i++)
            for (n = 0; n < *nw; n++)
                dN1[l*(*ns)*(*nw)+i*(*nw)+n] = 0.;

    for (n = 0; n < *nw; n++) {
        i=1-1; //N1
        l=0; li=l+i*2+n*2*(*ns);
        dN1[l*(*ns)*(*nw)+i*(*nw)+n] = -1.;
        l=1; li=l+i*2+n*2*(*ns);
        dN1[l*(*ns)*(*nw)+i*(*nw)+n] = -1.;

        i=2-1; //N2
        l=0; li=l+i*2+n*2*(*ns);
        dN1[l*(*ns)*(*nw)+i*(*nw)+n] = 1.;
        l=1; li=l+i*2+n*2*(*ns);
        dN1[l*(*ns)*(*nw)+i*(*nw)+n] = 0.;

        i=3-1; //N3
        l=0; li=l+i*2+n*2*(*ns);
        dN1[l*(*ns)*(*nw)+i*(*nw)+n] = 0.;
        l=1; li=l+i*2+n*2*(*ns);
        dN1[l*(*ns)*(*nw)+i*(*nw)+n] = 1.;
    }
}

/*poisson bracket tensor
  Eij[i][j][n] = dNi/dksi dbNj/deta - dNi/deta dbNj/dksi , i=1,18, j=1,3
*/
    for (i = 0; i < *nvar; i++)
    for (j = 0; j < *ns; j++)
        for (n = 0; n < *nw; n++) {
            Eij[i*(*ns)*(*nw)+j*(*nw)+n] = 
                           dN1[0*(*nvar)*(*nw)+i*(*nw)+n] * dbN1[1*(*ns)*(*nw) + j*(*nw) + n]
                         - dN1[1*(*nvar)*(*nw)+i*(*nw)+n] * dbN1[0*(*ns)*(*nw) + j*(*nw) + n];
    }

/*P bracket tensor
  Pijk[i][j][k] = NN[i] * Eij[j][k], rhs[i]=sum_j sum_k uu_j y_k P[i][j][k]
    for (i = 0; i < *nvar; i++)
       for (j = 0; j < *ns; j++)
       for (k = 0; k < *ns; k++) {
           Pijk[i*(*ns)*(*ns)+j*(*ns)+k]=0.;
           for (n = 0; n < *nw; n++)
           Pijk[i*(*ns)*(*ns)+j*(*ns)+k] += weight[n] * NN[i*(*nw)+n] * Eij[j*(*ns)*(*nw)+k*(*nw)+n];
       }
*/

    free_dVector1D(bNN , 0, (*ns)*(*nw) - 1);
    free_dVector1D(dbN1 , 0, (*ns)*(*ns)*(*nw) - 1);

  return 0;
}


int lagrange_element_matrix_(/*int *node,*/
                          int *ns,
                          int *nw,
                          int *ndof,
                          int *nvar,
                          double *x,
                          double *y,
                          double *weight,
                          /*double *dbN1,
                          double *dbN2,*/
                          double *NN,
                          double *dN1,
                          double *Eij,
                          double *bEij,
                          double *ma,
                          double *sa,
                          double *srr,
                          double *srz,
                          double *szr,
                          double *szz,
                          double *spp,
                          double *sr,
                          double *sz
                         )
{
/* matrix related, changes over each trianlge
*/
  int i, j, k, n, irow, ii, jj;
  double *jacobian;
  double **dnx, **dny;
  double **m, **s, **mtmp, **stmp, **tij, **tji;
  double **r, **z, **rtmp, **ztmp;
  double **rr, **rz, **rrtmp, **rztmp;
  double **zr, **zz, **zrtmp, **zztmp;
  double **pp, **pptmp;
  double J11, J12, J21, J22;
  double J31, J32, J41, J42, J51, J52;
  double ravg, area;
  static int count=0, pcount=0, ncount=0, parea=0, narea=0;


  jacobian = dVector1D(            0, *nw  -1);
  dnx      = dVector2D(0, *nvar-1, 0, *nw  -1);
  dny      = dVector2D(0, *nvar-1, 0, *nw  -1);
  m        = dVector2D(0, *nvar-1, 0, *nvar-1);
  s        = dVector2D(0, *nvar-1, 0, *nvar-1);
  mtmp     = dVector2D(0, *nvar-1, 0, *nvar-1);
  stmp     = dVector2D(0, *nvar-1, 0, *nvar-1);
  r        = dVector2D(0, *nvar-1, 0, *nvar-1);
  z        = dVector2D(0, *nvar-1, 0, *nvar-1);
  rtmp     = dVector2D(0, *nvar-1, 0, *nvar-1);
  ztmp     = dVector2D(0, *nvar-1, 0, *nvar-1);
  rr       = dVector2D(0, *nvar-1, 0, *nvar-1);
  rz       = dVector2D(0, *nvar-1, 0, *nvar-1);
  rrtmp    = dVector2D(0, *nvar-1, 0, *nvar-1);
  rztmp    = dVector2D(0, *nvar-1, 0, *nvar-1);
  zr       = dVector2D(0, *nvar-1, 0, *nvar-1);
  zz       = dVector2D(0, *nvar-1, 0, *nvar-1);
  zrtmp    = dVector2D(0, *nvar-1, 0, *nvar-1);
  zztmp    = dVector2D(0, *nvar-1, 0, *nvar-1);
  pp       = dVector2D(0, *nvar-1, 0, *nvar-1);
  pptmp    = dVector2D(0, *nvar-1, 0, *nvar-1);
  tij      = dVector2D(0, *nvar-1, 0, *nvar-1);
  tji      = dVector2D(0, *nvar-1, 0, *nvar-1);

/* jacobian over each triangle
*/
  for (n = 0; n < *nw; ++n) {
      jacobian[n] = 0.;
      for (i = 0; i < *ns; ++i)
      for (j = 0; j < *ns; ++j)
          jacobian[n] += bEij[i*(*ns)*(*nw)+j*(*nw)+n] * x[i] * y[j];
//      if(jacobian[n] < 0.) printf("jacobian[%d]=%e (<0) at e[%d]\n", n, jacobian[n], count);
  }
/*
  if(*ns>3) {
      if(jacobian[0]<0.||jacobian[1]<0.||jacobian[2]<0.||jacobian[3]<0.||jacobian[4]<0.||jacobian[5]<0.) ncount++;
      else pcount++;
      area=(x[2]-x[0])*(y[4]-y[2])-(x[4]-x[2])*(y[2]-y[0]);
      if(area<0.) narea++; else parea++;
      printf("jacobian<0 count=%d ; jacobian>0 count=%d area=%e\n", 
             ncount,pcount,narea,parea);
  } else {
      if(jacobian[0]<0.||jacobian[1]<0.||jacobian[2]<0.) ncount++;
      else pcount++;
      area=(x[1]-x[0])*(y[2]-y[1])-(x[2]-x[1])*(y[1]-y[0]);
      if(area<0.) narea++; else parea++;
      printf("jacobian<0 count=%d ; jacobian>0 count=%d area<0 count=%d; >0 count=%d\n", 
             ncount,pcount,narea,parea);
  }
*/

/* first order derivatives
 */
  for (n = 0; n < *nw; ++n)
      for (i = 0; i < *nvar; ++i) {
          dnx[i][n] = 0.;
          dny[i][n] = 0.;
          for (j = 0; j < *ns; ++j) {
              dnx[i][n] += Eij[i*(*ns)*(*nw)+j*(*nw)+n] * y[j] / jacobian[n];
              dny[i][n] -= Eij[i*(*ns)*(*nw)+j*(*nw)+n] * x[j] / jacobian[n];
          }
      }

/* average (r,theta) coordinate jacobian
 */
#ifdef RTHETA
  ravg = 0.;
  for (i = 0; i < *ns; ++i) ravg += x[i];
  ravg /= (*ns);
#else
  ravg = 1.;
#endif

/* transformation matrix
 */
//#ifdef ORDER2
//#else
  for (i = 0; i < *nvar; ++i)
      for (j = 0; j < *nvar; ++j) {
          tij[i][j]=0.;
          tji[i][j]=0.;
      }
  //at node i
  for (i = 0; i < *ns; ++i) {
      J11 = x[2-1] - x[1-1];
      J12 = y[2-1] - y[1-1];
      J21 = x[3-1] - x[1-1];
      J22 = y[3-1] - y[1-1];

      irow= i*(*ndof);
      tij[irow+0][irow+0]=1.;
      if(*ndof>1) {
      tij[irow+1][irow+1]=J11; tij[irow+1][irow+1+1]=J12;
      tij[irow+2][irow+1]=J21; tij[irow+2][irow+1+1]=J22;
      tij[irow+3][irow+3]=J11*J11; tij[irow+3][irow+3+1]=2.*J11*J12;      tij[irow+3][irow+3+2]=J12*J12; 
      tij[irow+4][irow+3]=J11*J21; tij[irow+4][irow+3+1]=J12*J21+J11*J22; tij[irow+4][irow+3+2]=J12*J22; 
      tij[irow+5][irow+3]=J21*J21; tij[irow+5][irow+3+1]=2.*J21*J22;      tij[irow+5][irow+3+2]=J22*J22; 
      }

      tji[irow+0][irow+0]=1.;
      if(*ndof>1) {
      tji[irow+1][irow+1]=J11; tji[irow+1+1][irow+1]=J12;
      tji[irow+1][irow+2]=J21; tji[irow+1+1][irow+2]=J22;
      tji[irow+3][irow+3]=J11*J11; tji[irow+4][irow+3+0]=2.*J11*J12;      tji[irow+5][irow+3+0]=J12*J12; 
      tji[irow+3][irow+4]=J11*J21; tji[irow+4][irow+3+1]=J12*J21+J11*J22; tji[irow+5][irow+3+1]=J12*J22; 
      tji[irow+3][irow+5]=J21*J21; tji[irow+4][irow+3+2]=2.*J21*J22;      tji[irow+5][irow+3+2]=J22*J22; 
      }
  }
//#endif

/* calculate mass (m), stiffness (s) matrices 
   over collocation points
*/
  for (j = 0; j < *nvar; ++j)
  for (i = 0; i < *nvar; ++i) {
          m[j][i]  = 0.;
          s[j][i]  = 0.;
          rr[j][i]  = 0.;
          rz[j][i]  = 0.;
          zr[j][i]  = 0.;
          zz[j][i]  = 0.;
          pp[j][i]  = 0.;
          r[j][i]  = 0.;
          z[j][i]  = 0.;
  for (n = 0; n < *nw; ++n) {
/*
#ifdef ORDER2
      for (ii = 0; ii < *nvar; ++ii)
      for (jj = 0; jj < *nvar; ++jj) {
          tij[ii][jj]=0.; tji[ii][jj]=0.;
      }
      J11=0.; J12=0.;
      J21=0.; J22=0.;
      J31=0.; J32=0.;
      J41=0.; J42=0.;
      J51=0.; J52=0.;
      for (j = 0; j < *ns; j++) {
          J11 += x[j] * dbN1[0*(*ns)*(*nw) + j*(*nw) + n];
          J12 += y[j] * dbN1[0*(*ns)*(*nw) + j*(*nw) + n];
          J21 += x[j] * dbN1[1*(*ns)*(*nw) + j*(*nw) + n];
          J22 += y[j] * dbN1[1*(*ns)*(*nw) + j*(*nw) + n];

          J31 += x[j] * dbN2[0*(*ns)*(*nw) + j*(*nw) + n];
          J32 += y[j] * dbN2[0*(*ns)*(*nw) + j*(*nw) + n];
          J41 += x[j] * dbN2[1*(*ns)*(*nw) + j*(*nw) + n];
          J42 += y[j] * dbN2[1*(*ns)*(*nw) + j*(*nw) + n];
          J51 += x[j] * dbN2[2*(*ns)*(*nw) + j*(*nw) + n];
          J52 += y[j] * dbN2[2*(*ns)*(*nw) + j*(*nw) + n];
      }
      for (k = 0; k < *node; ++k) {
          irow= k*(*ndof);
          tij[irow+0][irow+0]=1.;
          tij[irow+1][irow+1]=J11; tij[irow+1][irow+1+1]=J12;
          tij[irow+2][irow+1]=J21; tij[irow+2][irow+1+1]=J22;
          tij[irow+3][irow+1]=J31; tij[irow+3][irow+1+1]=J32;
          tij[irow+4][irow+1]=J41; tij[irow+4][irow+1+1]=J42;
          tij[irow+5][irow+1]=J51; tij[irow+5][irow+1+1]=J52;
          tij[irow+3][irow+3]=J11*J11; tij[irow+3][irow+3+1]=2.*J11*J12;      tij[irow+3][irow+3+2]=J12*J12; 
          tij[irow+4][irow+3]=J11*J21; tij[irow+4][irow+3+1]=J12*J21+J11*J22; tij[irow+4][irow+3+2]=J12*J22; 
          tij[irow+5][irow+3]=J21*J21; tij[irow+5][irow+3+1]=2.*J21*J22;      tij[irow+5][irow+3+2]=J22*J22; 

          tji[irow+0][irow+0]=1.;
          tji[irow+1][irow+1]=J11; tji[irow+1+1][irow+1]=J12;
          tji[irow+1][irow+2]=J21; tji[irow+1+1][irow+2]=J22;
          tji[irow+1][irow+3]=J31; tji[irow+1+1][irow+3]=J32;
          tji[irow+1][irow+4]=J41; tji[irow+1+1][irow+4]=J42;
          tji[irow+1][irow+5]=J51; tji[irow+1+1][irow+5]=J52;
          tji[irow+3][irow+3]=J11*J11; tji[irow+4][irow+3+0]=2.*J11*J12;      tji[irow+5][irow+3+0]=J12*J12; 
          tji[irow+3][irow+4]=J11*J21; tji[irow+4][irow+3+1]=J12*J21+J11*J22; tji[irow+5][irow+3+1]=J12*J22; 
          tji[irow+3][irow+5]=J21*J21; tji[irow+4][irow+3+2]=2.*J21*J22;      tji[irow+5][irow+3+2]=J22*J22; 
      }
      for (jj = 0; jj < *nvar; ++jj)
      for (ii = 0; ii < *nvar; ++ii) {
          m[  j][i] += tij[jj][i]*NN[jj*(*nw)+n] * NN[ii*(*nw)+n]*tij[ii][i];
          m[  j][i] *= weight[n] * jacobian[n];
          s[  j][i] += tij[jj][i]*(dnx[jj][n]*dnx[ii][n] +  dny[jj][n]*dny[ii][n])*tij[ii][i];
          s[  j][i] *= weight[n] * jacobian[n];
      }
#else
#endif
*/
#ifdef RTHETA
          m[  j][i] += weight[n] * jacobian[n] * NN[j*(*nw)+n] * NN[i*(*nw)+n]*ravg;
          s[  j][i] -= weight[n] * jacobian[n] * (dnx[j][n]*dnx[i][n]*ravg + dny[j][n]*dny[i][n]/ravg);
          rr[j][i] -= weight[n] * jacobian[n] * (dnx[j][n]*dnx[i][n]*ravg);
          rz[j][i] -= weight[n] * jacobian[n] * (dnx[j][n]*dny[i][n]);
          zr[j][i] -= weight[n] * jacobian[n] * (dny[j][n]*dnx[i][n]);
          zz[j][i] -= weight[n] * jacobian[n] * (dny[j][n]*dny[i][n]/ravg);
          pp[j][i] -= weight[n] * jacobian[n] * (dny[j][n]*dny[i][n]);
          r[  j][i] += weight[n] * jacobian[n] * (NN[i*(*nw)+n]*dnx[j][n]*ravg);
          z[  j][i] += weight[n] * jacobian[n] * (NN[i*(*nw)+n]*dny[j][n]/ravg);
#else
          m[  j][i] += weight[n] * jacobian[n] * NN[j*(*nw)+n] * NN[i*(*nw)+n];
          s[  j][i] -= weight[n] * jacobian[n] * (dnx[j][n]*dnx[i][n] + dny[j][n]*dny[i][n]);
          rr[j][i] -= weight[n] * jacobian[n] * (dnx[j][n]*dnx[i][n]);
          rz[j][i] -= weight[n] * jacobian[n] * (dnx[j][n]*dny[i][n]);
          zr[j][i] -= weight[n] * jacobian[n] * (dny[j][n]*dnx[i][n]);
          zz[j][i] -= weight[n] * jacobian[n] * (dny[j][n]*dny[i][n]);
          pp[j][i] -= weight[n] * jacobian[n] * (dny[j][n]*dny[i][n]);
          r[ j][i] += weight[n] * jacobian[n] * (NN[i*(*nw)+n]*dnx[j][n]);
          z[ j][i] += weight[n] * jacobian[n] * (NN[i*(*nw)+n]*dny[j][n]);
#endif
      }
//  PetscPrintf(PETSC_COMM_WORLD,
//              "\t j=%d i=%d m=%e s=%e\n", j, i, m[j][i], s[j][i]);
  }

//#ifdef ORDER2
//#else
  for (j = 0; j < *nvar; ++j)
  for (i = 0; i < *nvar; ++i) {
          mtmp[j][i] = 0.;
          stmp[j][i] = 0.;
          rrtmp[j][i] = 0.;
          rztmp[j][i] = 0.;
          zrtmp[j][i] = 0.;
          zztmp[j][i] = 0.;
          pptmp[j][i] = 0.;
          rtmp[j][i] = 0.;
          ztmp[j][i] = 0.;
      for (k = 0; k < *nvar; ++k) {
          mtmp[j][i] += m[j][k] * tij[k][i];
          stmp[j][i] += s[j][k] * tij[k][i];
          rrtmp[j][i] += rr[j][k] * tij[k][i];
          rztmp[j][i] += rz[j][k] * tij[k][i];
          zrtmp[j][i] += zr[j][k] * tij[k][i];
          zztmp[j][i] += zz[j][k] * tij[k][i];
          pptmp[j][i] += pp[j][k] * tij[k][i];
          rtmp[j][i] += r[j][k] * tij[k][i];
          ztmp[j][i] += z[j][k] * tij[k][i];
      }
  }

  for (j = 0; j < *nvar; ++j)
  for (i = 0; i < *nvar; ++i) {
          ma[i+j*(*nvar)] = 0.;
          sa[i+j*(*nvar)] = 0.;
          srr[i+j*(*nvar)] = 0.;
          srz[i+j*(*nvar)] = 0.;
          szr[i+j*(*nvar)] = 0.;
          szz[i+j*(*nvar)] = 0.;
          spp[i+j*(*nvar)] = 0.;
          sr[i+j*(*nvar)] = 0.;
          sz[i+j*(*nvar)] = 0.;
      for (k = 0; k < *nvar; ++k) {
          ma[i+j*(*nvar)] += tji[j][k] * mtmp[k][i];
          sa[i+j*(*nvar)] += tji[j][k] * stmp[k][i];
          srr[i+j*(*nvar)] += tji[j][k] * rrtmp[k][i];
          srz[i+j*(*nvar)] += tji[j][k] * rztmp[k][i];
          szr[i+j*(*nvar)] += tji[j][k] * zrtmp[k][i];
          szz[i+j*(*nvar)] += tji[j][k] * zztmp[k][i];
          spp[i+j*(*nvar)] += tji[j][k] * pptmp[k][i];
          sr[i+j*(*nvar)] += tji[j][k] * rtmp[k][i];
          sz[i+j*(*nvar)] += tji[j][k] * ztmp[k][i];
      }
//  PetscPrintf(PETSC_COMM_WORLD,
//              "\t j=%d i=%d m=%e s=%e\n", j, i, ma[i+j*(*nvar)], sa[i+j*(*nvar)]);
  }
//#endif

  free_dVector1D(jacobian,         0, *nw-1);
  free_dVector2D(dnx,  0, *nvar-1, 0, *nw-1);
  free_dVector2D(dny,  0, *nvar-1, 0, *nw-1);
  free_dVector2D(m,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(s,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(mtmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(stmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(r,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(z,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(rtmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(ztmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(rr,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(rz,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(rrtmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(rztmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(zr,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(zz,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(zrtmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(zztmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(pp,    0, *nvar-1, 0, *nvar-1);
  free_dVector2D(pptmp, 0, *nvar-1, 0, *nvar-1);
  free_dVector2D(tij,  0, *nvar-1, 0, *nvar-1);
  free_dVector2D(tji,  0, *nvar-1, 0, *nvar-1);
  if(count==0) PetscPrintf(PETSC_COMM_WORLD,
              "\ttric1_element_matrix: done\n");
  count++;
  return 0;
}
