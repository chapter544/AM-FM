/*
 *  Function: recon2d_c
 *
 *  Function returns the reconstructed real and imaginary image components
 *    computed from the AM and FM components and the 4 initial phase values.
 *
 *  10/21/2005  ras
 *
 */

/* The matlab implementation of this is actually nearly as fast as the 
 * C implementation */

#include <mex.h>
#include <math.h>

#define PHASE_SCALE 300

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* function declarations */
short IsplineGrad(double *dx,double *dy,double *img,int numRows,int numCols,
                  short horizSym,short vertSym,
                  double imgUL,double imgUR,double imgLL,double imgLR);
short IsplineDeriv(double *xdot,double *x,int len,short sym,double x0,
                   double x1);





/* ----------------------------------------------------------------------------*/
/*  
 *  Function: mexFunction
 *
 *  equivilant to matlab function 
 *  function [img, p] = recon2d(a, u, v, pInit, numRows, numCols)
 */
void mexFunction(int nlhs,	     /* Num return vals on lhs */
		 mxArray *plhs[],    /* Matrices on lhs      */
		 int nrhs,	     /* Num args on rhs    */
		 const mxArray *prhs[]     /* Matrices on rhs */
		 )
  {
  double *a, *u, *v, *pInit, *imgR, *imgI, *p;
  int numRows, numCols;
  int            pixCount;              /* pixel counter */
  int            numPix;                /* number of pixels in the image */


  if (nrhs != 6) mexErrMsgTxt("Requires arguments (a, u, v, pInit, numRows, numCols).");

  /*read parameters */
  a = mxGetPr(prhs[0]);
  u = mxGetPr(prhs[1]);
  v = mxGetPr(prhs[2]);
  pInit = mxGetPr(prhs[3]);
  numRows = (int)(mxGetScalar(prhs[4]));
  numCols = (int)(mxGetScalar(prhs[5]));

  /*allocate output arrays */
  plhs[0] = (mxArray *) mxCreateDoubleMatrix(numCols,numRows,mxCOMPLEX);
  if (plhs[0] == NULL) mexErrMsgTxt("Error allocating result matrix");
  imgR = mxGetPr(plhs[0]);
  imgI = mxGetPi(plhs[0]);
  plhs[1] = (mxArray *) mxCreateDoubleMatrix(numCols,numRows,mxREAL);
  if (plhs[1] == NULL) mexErrMsgTxt("Error allocating result matrix");
  p = mxGetPr(plhs[1]);



 
 
  numPix = numRows*numCols;

  /* allocate the working arrays 
  if((p = (double *)malloc(numPix*sizeof(double))) == NULL) {
    printf("\nrecon2d: Free store exhausted.\n");
    return -1;
  }*/

  /* scale the phase gradient */
  for(pixCount = 0; pixCount < numPix; pixCount++) {
    u[pixCount] *= PHASE_SCALE;
    v[pixCount] *= PHASE_SCALE;
  }

  /* reconstruct the phase modulation */
  /* Note: seeing same thing as in demod where u and v need to switch 
   * here that means switching directions in pInit as well 
   * haven't nailed this down yet */
  if(IsplineGrad(v,u,p,numRows,numCols,1,1,
                 pInit[0],pInit[2],pInit[1],pInit[3]) != 0) {
    printf("\nrecon2d: Phase reconstruction failed.\n");
    free(p);
    return;
  }

  /* reconstruct the real and imaginary components */
  for(pixCount = 0; pixCount < numPix; pixCount++) {
    imgR[pixCount] = a[pixCount]*cos(p[pixCount]);
    imgI[pixCount] = a[pixCount]*sin(p[pixCount]);
  }


  return;
  }      
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: IsplineGrad
 *
 *  Function computes the original image img that was used by splineGrad
 *    to produce the gradient images dx and dy. numRows and numCols are
 *    the number of image rows and columns, respectively. If horizontal or
 *    vertical mirror symmetry or no symmetry boundary conditions are imposed,
 *    then the four pixels from the upper left corner of the image are required
 *    initial conditions and are given in imgUL, imgUR, imgLL, and imgLR.
 *    Returns -1 on error.
 *
 *  5/28/2004  ras
 *
 *  REVISED 4/25/2005 ras: added error messages for reconstruction failure.
 *
 *  REVISED 8/12/2005 ras: added horizontal and vertical mirror antisymmetry and
 *    no symmetry boundary conditions for the output data.
 *
 *    Symmetry		horizSym/vertSym
 *    --------		----------------
 *    None		0
 *    Mirror		positive
 *    Mirror Anti	negative
 *
 */

short IsplineGrad(double *dx,double *dy,double *img,int numRows,int numCols,
                  short horizSym,short vertSym,
                  double imgUL,double imgUR,double imgLL,double imgLR) {
  double       *colData;            /* temp column data array */
  double       *gradData;           /* temp gradient data array */
  int           row;                /* row counter */
  int           col;                /* col counter */
  
  /* check input size */
  if(numRows < 2) {
    printf("\nIsplineGrad: Insufficient number of rows.\n");
    return -1;
  }
  if(numCols < 2) {
    printf("\nIsplineGrad: Insufficient number of columns.\n");
    return -1;
  }
  
  /* allocate temporary column and gradient arrays */
  if((colData = (double *)malloc(numRows*sizeof(double))) == NULL) {
    printf("\nIsplineGrad: Free store exhausted.\n");
    return -1;
  }
  if((gradData = (double *)malloc(numRows*sizeof(double))) == NULL) {
    printf("\nIsplineGrad: Free store exhausted.\n");
    free(colData);
    return -1;
  }
  
  /* reconstruct the first two columns from the vertical gradient */
  col = 0;
  for(row = 0; row < numRows; row++) {
    gradData[row] = dy[row*numCols + col];
  }
  if(IsplineDeriv(gradData,colData,numRows,vertSym,imgUL,imgLL) != 0) {
    printf("\nIsplineGrad: Reconstruction failed on column %d.\n",col);
    free(colData);
    free(gradData);
    return -1;
  }
  for(row = 0; row < numRows; row++) {
    img[row*numCols + col] = colData[row];
  }
  col = 1;
  for(row = 0; row < numRows; row++) {
    gradData[row] = dy[row*numCols + col];
  }
  if(IsplineDeriv(gradData,colData,numRows,vertSym,imgUR,imgLR) != 0) {
    printf("\nIsplineGrad: Reconstruction failed on column %d.\n",col);
    free(colData);
    free(gradData);
    return -1;
  }
  for(row = 0; row < numRows; row++) {
    img[row*numCols + col] = colData[row];
  }
  
  free(colData);
  free(gradData);
  
  /* reconstruct the rows from the horizontal gradient */
  for(row = 0; row < numRows; row++) {
    if(IsplineDeriv(dx+row*numCols,img+row*numCols,numCols,horizSym,
                    img[row*numCols],img[row*numCols + 1]) != 0) {
      printf("\nIsplineGrad: Reconstruction failed on row %d.\n",row);
      return -1;
    }
  }
  
  return 0;
}
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: IsplineDeriv
 *
 *  Function computes the original values x that were used by splineDeriv
 *    to produce xdot. If no symmetry or mirror-symmetric boundary conditions
 *    are imposed, then x0 and x1 are required initial conditions for perfect
 *    reconstruction. See
 *      (1) "B-Spline Signal Processing: Part I..." pp. 824,826-827
 *      (2) "B-Spline Signal Processing: Part II..." pp. 835
 *    Returns -1 on error.
 *
 *  5/28/2004  ras
 *
 *  REVISED 8/12/2005 ras: added mirror antisymmetry and no symmetry boundary
 *    conditions for the output data. See pages 117 through 123 of thesis
 *    notebook.
 *
 *    Symmetry		sym
 *    --------		---
 *    None		0
 *    Mirror		positive
 *    Mirror Anti	negative
 *
 */

short IsplineDeriv(double *xdot,double *x,int len,short sym,double x0,
                   double x1) {
  int           i;                  /* loop counter */
  
  /* check input data length */
  if(len < 2) {
    printf("\nIsplineDeriv: Insufficient data length.\n");
    return -1;
  }
  
  /* reconstruct the original samples from the derivative */
  if(sym < 0) {                   /* mirror antisymmetry */

    x[0] = 0;
    x[1] = (xdot[1] + 2*xdot[0])/3;
  } else {                        /* mirror symmetry and no symmetry */
    x[0] = x0;
    x[1] = x1;
  }
  for(i = 2; i < len; i++) {
    x[i] = (xdot[i] + 4*xdot[i-1] + xdot[i-2])/3 + x[i-2];
  }
  
  return 0;
}
/* ----------------------------------------------------------------------------*/





