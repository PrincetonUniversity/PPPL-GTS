/*****************************************************************************/
/*                                                                           */
/*  (mesh_gen.c)                                                             */
/*                                                                           */
/*****************************************************************************/

/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#error
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <stdio.h>
#include <stdlib.h>
#include "triangle.h"
/*#include <petsc.h> */

/*****************************************************************************/
/*                                                                           */
/* mesh_gen () Create triangulation                                          */
/*                                                                           */
/*  Input: (graph)                                                           */
/*   npoints:                                                                */
/*   points:                                                                 */
/*   nedges:                                                                 */
/*   edges:                                                                  */
/*   flag: 0 for create and strore mesh (return 'ntri'); 1 for return mesh   */
/*                                                                           */
/*  Output:                                                                  */
/*   ntri: number of trianles made                                           */
/*   tri1[ntri], tri2[ntri], tri3[ntri]                                      */
/*                                                                           */
/*****************************************************************************/

/* PETSC_EXTERN_CXX_BEGIN */

static int *tries = NULL;
//static REAL *dpoints = NULL;

//void mesh_gen_( int *ntri, int *tri1, int *tri2, int *tri3 )
void mesh_gen_( int *order, int *ntri, int *tri)
{
  int i,j; 
  /* output */
  for( i = 0, j = 0; i < *ntri; i++ ) {
//    tri1[i] = tries[j++]; tri2[i] = tries[j++]; tri3[i] = tries[j++];
    if(*order == 2) {
    tri[i*6+0] = tries[j++];
    tri[i*6+2] = tries[j++];
    tri[i*6+4] = tries[j++];
    tri[i*6+3] = tries[j++];
    tri[i*6+5] = tries[j++];
    tri[i*6+1] = tries[j++];
    } else {
    tri[i*3+0] = tries[j++];
    tri[i*3+1] = tries[j++];
    tri[i*3+2] = tries[j++];
    }
  }
  free( tries );
  tries = NULL;
}

void mesh_count_( int *npoints, REAL *points, int *nedges, int *edges, int *order, int *ngvert, int *ntri, int *mype )
{
  struct triangulateio in, out, vorout;
  REAL hole[2] = {0., 0.};
  int i,j;
  
  /* Define input points. */
  in.numberofpoints = *npoints;
  in.numberofpointattributes = 0;
  in.pointlist = points;
  in.pointmarkerlist = (int *) NULL;   
  in.pointattributelist = (REAL *) NULL;   
    
  in.trianglelist = (int *) NULL;                          /* In / out */
  in.triangleattributelist = (REAL *) NULL;                /* In / out */
  in.trianglearealist = (REAL *) NULL;                     /* In only */
  
  in.numberoftriangles = 0;                                /* In / out */
  in.numberofcorners = 0;                                  /* In / out */
  in.numberoftriangleattributes - 0;                       /* In / out */
  
  in.regionlist = (REAL *) NULL;   /* In / pointer to array copied out */
  
  in.numberofholes = 1;
  in.numberofregions = 0;
  in.holelist = hole;
  
  in.segmentlist = edges;
  in.numberofsegments = *nedges;
  in.segmentmarkerlist = (int *) NULL;
  
  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out' and a voronoi diagram in `vorout'.  */
  out.trianglelist = (int *) NULL;
  out.numberoftriangleattributes = 0;
  out.neighborlist = (int *) NULL;     /* Needed only if -n switch used. */
  
  /* Triangulate the points.  Switches are chosen to read and write a  */
  /*   PSLG (p), ignore BCs (B), number everything from  */
  /*   zero (z), no node output (N), and verbose/quit (V/Q) */
  if(*order == 2) {
  triangulate( "pPNBQo2", &in, &out, &vorout );
  } else {
  triangulate( "pPNBQ", &in, &out, &vorout );
  }

if(*mype==0){
  /* file output*/
  FILE *outfile = fopen("n.node", "w");
  fprintf(outfile, "%ld  %d  %d  %d\n",in.numberofpoints,2,0,0);
  for(i=0;i<in.numberofpoints * 2;i += 2){
    fprintf(outfile, "%4d    %.17g  %.17g\n",i/2 + 1, points[i],points[i+1]);
  }
  fclose(outfile);
  outfile = fopen("n.ele", "w");
     if(*order == 2) {
  fprintf(outfile, "%ld  %d  %d\n", out.numberoftriangles, 6, 0);
  for(j=0;j<out.numberoftriangles*6;j += 6){
    fprintf(outfile, "%4ld    %4d  %4d  %4d  %4d  %4d  %4d\n",j/6+1,
            out.trianglelist[j+0], out.trianglelist[j+1], out.trianglelist[j+2],
            out.trianglelist[j+3], out.trianglelist[j+4], out.trianglelist[j+5]
           );
  }
     } else {
  fprintf(outfile, "%ld  %d  %d\n", out.numberoftriangles, 3, 0);
  for(j=0;j<out.numberoftriangles*3;j += 3){
    fprintf(outfile, "%4ld    %4d  %4d  %4d\n",j/3+1,
            out.trianglelist[j], out.trianglelist[j+1], out.trianglelist[j+2]);
  }
     }
  fclose(outfile);
};
    
  /* output */
  *ntri = out.numberoftriangles;
  *ngvert = out.numberofpoints;
  /* save */
  tries = out.trianglelist;
//  if(*order > 1) dpoints = out.pointlist;
 }

/* PETSC_EXTERN_CXX_END */
