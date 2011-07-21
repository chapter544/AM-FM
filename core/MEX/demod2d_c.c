/*
 *  Function: demod2d_c
 *
 *  Function demodulates a complex image into its AM and FM components. The FM
 *    component is computed using a tensor product cubic spline model applied to
 *    the scaled, unwrapped phase embedded with congruent principal values. This
 *    unwrapped phase modulation is also returned. At pixels where the amplitude
 *    is zero, the phase is copied from a neighboring pixel with a valid
 *    principal value.
 *
 *  10/21/2005  ras
 *
 */

#include <fftw3.h>
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
short unwrap2d(double *ppv,double *p,int numRows,int numCols);
double xround(double x);
short phaseGrad2est(double *ppv,double *u2,double *v2,int numRows,int numCols,
                    short horizSym,short vertSym);
short splineGrad2(double *img,double *dx2,double *dy2,int numRows,int numCols,
                  short horizSym,short vertSym);
short splineDeriv2(double *x,double *xdotdot,int len,short sym);
short splineGrad(double *img,double *dx,double *dy,int numRows,int numCols,
                 short horizSym,short vertSym);
short splineDeriv(double *x,double *xdot,int len,short sym);
short symExpFilt(double *x,double *cMinus,int len,short sym,double pole,
                 double gain);





/* ----------------------------------------------------------------------------*/
/*  
 *  Function: mexFunction
 *
 *  equivilant to matlab function 
 *  function [a, u, v, p] = demod2d_c(img, numRows, numCols)
 */
void mexFunction(int nlhs,	     /* Num return vals on lhs */
		 mxArray *plhs[],    /* Matrices on lhs      */
		 int nrhs,	     /* Num args on rhs    */
		 const mxArray *prhs[]     /* Matrices on rhs */
		 )
  {
  double *imgR, *imgI, *a, *u, *v, *p;
  int numRows, numCols;
  int            row;                   /* row counter; */
  int            col;                   /* col counter; */
  int            pixCount;              /* pixel counter */
  int            numPix;                /* number of pixels in the image */
  int            numZeroAM;             /* number of pixels with zero amplitude */
  short         *zeroAM;                /* 1 for zero amplitude, 0 otherwise */

  if (nrhs != 3) mexErrMsgTxt("Requires arguments (img,numRows,numCols).");

  /*read parameters */
  imgR = mxGetPr(prhs[0]);
  imgI = mxGetPi(prhs[0]);
  numRows = (int)(mxGetScalar(prhs[1]));
  numCols = (int)(mxGetScalar(prhs[2]));

  /*allocate output arrays */
  plhs[0] = (mxArray *) mxCreateDoubleMatrix(numCols,numRows,mxREAL);
  if (plhs[0] == NULL) mexErrMsgTxt("Error allocating result matrix");
  a = mxGetPr(plhs[0]);

  plhs[1] = (mxArray *) mxCreateDoubleMatrix(numCols,numRows,mxREAL);
  if (plhs[1] == NULL) mexErrMsgTxt("Error allocating result matrix");
  v = mxGetPr(plhs[1]);

  plhs[2] = (mxArray *) mxCreateDoubleMatrix(numCols,numRows,mxREAL);
  if (plhs[2] == NULL) mexErrMsgTxt("Error allocating result matrix");
  u = mxGetPr(plhs[2]);

  /* Note: One would think, according to the function definitions, that 
     u would go to pointer 1 and v to pointer 2, but doing that causes 
     the output images to be switched and I've not yet figured out why */
  plhs[3] = (mxArray *) mxCreateDoubleMatrix(numCols,numRows,mxREAL);
  if (plhs[3] == NULL) mexErrMsgTxt("Error allocating result matrix");
  p = mxGetPr(plhs[3]);


  numPix = numRows*numCols;

  /* allocate the working arrays */
  if((zeroAM = (short *)malloc(numPix*sizeof(short))) == NULL) {
    printf("\ndemod2d: Free store exhausted.\n");
    return;
  }

  /* compute the amplitude and phase modulations */
  numZeroAM = 0;
  for(pixCount = 0; pixCount < numPix; pixCount++) {
    a[pixCount] = sqrt(imgR[pixCount]*imgR[pixCount] +
                       imgI[pixCount]*imgI[pixCount]);
    if((fabs(imgR[pixCount]) < DBL_EPSILON) &&
       (fabs(imgI[pixCount]) < DBL_EPSILON)) {
      /* set zeroAM flag and increment its counter */
      zeroAM[pixCount] = 1;
      numZeroAM++;
    } else {
      /* clear zeroAM flag */
      zeroAM[pixCount] = 0;
      /* calculate wrapped phase in (-pi,pi] */
      p[pixCount] = atan2(imgI[pixCount],imgR[pixCount]);
      if(p[pixCount] == -M_PI) {
        p[pixCount] = M_PI;
      }
    }
  }

  /* return when amplitude is zero everywhere */
  if(numZeroAM == numPix) {
    printf("\ndemod2d: Cannot demodulate image with no support.\n");
    free(zeroAM);
    return;
  }

  /* interpolate the phase at pixels with zero amplitude by copying a valid */
  /*   phase value from a neighboring pixel; neighbors are searched in the order */
  /*   left, right, top, bottom, and mirror symmetry is used at the boundaries */
  
  while(numZeroAM > 0) {
    for(row = 0; row < numRows; row++) {
      for(col = 0; col < numCols; col++) {
        if(zeroAM[row*numCols + col] != 0) {
          /* process the left phase value */
          if(col == 0) {
            if(zeroAM[row*numCols + col+1] == 0) {
              p[row*numCols + col] = p[row*numCols + col+1];
              zeroAM[row*numCols + col] = 0;
              numZeroAM--;
              continue;
            }
          } else {
            if(zeroAM[row*numCols + col-1] == 0) {
              p[row*numCols + col] = p[row*numCols + col-1];
              zeroAM[row*numCols + col] = 0;
              numZeroAM--;
              continue;
            }
          }
          /* process the right phase value */
          if(col == numCols-1) {
            if(zeroAM[row*numCols + col-1] == 0) {
              p[row*numCols + col] = p[row*numCols + col-1];
              zeroAM[row*numCols + col] = 0;
              numZeroAM--;
              continue;
            }
          } else {
            if(zeroAM[row*numCols + col+1] == 0) {
              p[row*numCols + col] = p[row*numCols + col+1];
              zeroAM[row*numCols + col] = 0;
              numZeroAM--;
              continue;
            }
          }
          /* process the top phase value */
          if(row == 0) {
            if(zeroAM[(row+1)*numCols + col] == 0) {
              p[row*numCols + col] = p[(row+1)*numCols + col];
              zeroAM[row*numCols + col] = 0;
              numZeroAM--;
              continue;
            }
          } else {
            if(zeroAM[(row-1)*numCols + col] == 0) {
              p[row*numCols + col] = p[(row-1)*numCols + col];
              zeroAM[row*numCols + col] = 0;
              numZeroAM--;
              continue;
            }
          }
          /* process the bottom phase value */
          if(row == numRows-1) {
            if(zeroAM[(row-1)*numCols + col] == 0) {
              p[row*numCols + col] = p[(row-1)*numCols + col];
              zeroAM[row*numCols + col] = 0;
              numZeroAM--;
              continue;
            }
          } else {
            if(zeroAM[(row+1)*numCols + col] == 0) {
              p[row*numCols + col] = p[(row+1)*numCols + col];
              zeroAM[row*numCols + col] = 0;
              numZeroAM--;
              continue;
            }
          }
        }
      }
    }
  }
  free(zeroAM);


  /* unwrap the phase modulation */
  if(unwrap2d(p,p,numRows,numCols) != 0) {
    printf("\ndemod2d: Phase unwrapping failed.\n");
    return;
  }

  /* compute the phase gradient using a tensor product cubic spline model */
  if(splineGrad(p,u,v,numRows,numCols,1,1) != 0) {
    printf("\ndemod2d: Frequency estimation failed.\n");
    return;
  }

  /* remove the phase scale from the gradient */
  for(pixCount = 0; pixCount < numPix; pixCount++) {
    u[pixCount] /= PHASE_SCALE;
    v[pixCount] /= PHASE_SCALE;
  }

  


  return;
  }      
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: unwrap2d
 *
 *  Function computes the unwrapped phase p that minimizes the least-squares
 *    error between its second (cubic spline) derivatives and the estimated
 *    second phase derivatives. The principal phase values are embedded in a
 *    scaled version of the LS solution to simultaneously preserve the phase
 *    congruency and the (scaled) phase derivatives. Mirror symmetry is used at
 *    the boundaries. Note that the pointers ppv and p can be the same memory
 *    address. Returns -1 on error.
 *
 *  10/21/2005  ras
 *
 */

short unwrap2d(double *ppv,double *p,int numRows,int numCols) {
  double        *data1;                 /* FFT/IFFT array 1 */
  double        *data2;                 /* FFT/IFFT array 2 */
  double         ARfilt;                /* auto-regressive filter */
  double         horizWeight;           /* horizontal filtering weight */
  double         vertWeight;            /* vertical filtering weight */
  int            row;                   /* row counter */
  int            col;                   /* col counter */
  int            numPix;                /* number of pixels */
  int            pixCount;              /* pixel counter */
  fftw_plan      plan1;                 /* FFTW plan for data1 => DCT(data1) */
  fftw_plan      plan2;                 /* FFTW plan for data2 => DCT(data2) */
 
  numPix = numRows*numCols;

  /* allocate working arrays */
  if((data1 = (double *)malloc(numPix*sizeof(double))) == NULL) {
    printf("\nunwrap2d: Free store exhausted.\n");
    return -1;
  }
  if((data2 = (double *)malloc(numPix*sizeof(double))) == NULL) {
    printf("\nunwrap2d: Free store exhausted.\n");
    free(data1);
    return -1;
  }

  /* build FFTW plans for DCT-I transforms of real data data1, data2 */
  /* i.e. Discrete Cosine Transform Type I, ( Logical N=2*(n-1) ) */
  /* note that this DCT is its own inverse (ignoring the normalization) */
  plan1 = fftw_plan_r2r_2d(numRows,numCols,data1,data1,
                           FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE);
  plan2 = fftw_plan_r2r_2d(numRows,numCols,data2,data2,
                           FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE);

  /* estimate the second phase derivatives using mirror symmetric boundary */
  /*   conditions */
  if(phaseGrad2est(ppv,data1,data2,numRows,numCols,1,1) != 0) {
    printf("\nunwrap2d: Second phase derivative estimation failed.\n");
    free(data1);
    free(data2);
    fftw_destroy_plan(plan1);
    fftw_destroy_plan(plan2);
    return -1;
  }

  /* compute the (real and mirror symmetric) spectrums of the second phase */
  /*   derivatives */
  fftw_execute(plan1);
  fftw_execute(plan2);

  /* filter to compute the (real, zero-mean, and mirror symmetric) least-squares */
  /*   phase spectrum; note that the filters are purely real and periodically */
  /*   even */
  for(row = 0; row < numRows; row++) {
    vertWeight = 6 - 18/(cos(2*M_PI*(double)row/(2*numRows-2)) + 2);
    for(col = 0; col < numCols; col++) {
      horizWeight = 6 - 18/(cos(2*M_PI*(double)col/(2*numCols-2)) + 2);
      ARfilt = horizWeight*horizWeight + vertWeight*vertWeight;
      if(fabs(ARfilt) < DBL_EPSILON) {
        /* the AR filter will be zero at dc; assume phase is zero mean */
        data1[row*numCols + col] = 0;
      } else {
        /* compute the spectrum of the phase */
        data1[row*numCols + col] = data1[row*numCols + col]*horizWeight +
                                   data2[row*numCols + col]*vertWeight;
        data1[row*numCols + col] /= ARfilt;
        /* normalize the transform */
        data1[row*numCols + col] /= 4*(numRows-1)*(numCols-1);
      }
    }
  }

  /* compute the LS phase from its spectrum */
  fftw_execute(plan1);
  
  /* embed the principal phase values in the scaled LS solution */
  for(pixCount = 0; pixCount < numPix; pixCount++) {
    p[pixCount] = ppv[pixCount] +
      2*M_PI*xround((data1[pixCount]*PHASE_SCALE - ppv[pixCount])/(2*M_PI));
  }

  free(data1);
  free(data2);
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  return 0;
}
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: xround
 *
 *  Rounds to the nearest integer
 *
 *  round() is not in standard math.h 
 *
 */
double xround(double x) {
  if(x < 0){
    return ceil(x-0.5);
  }else{
    return floor(x+0.5);
  }
}
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: phaseGrad2est
 *
 *  Function returns the second horizontal and vertical phase derivatives u2,v2
 *    estimated from the principal value phase ppv. Derivatives are computed
 *    using cubic spline models with boundary conditions given by horizSym
 *    and vertSym. Returns -1 on error.
 *
 *  5/17/2005  ras
 *
 *  REVISED 8/15/2005 ras: added horizontal and vertical symmetry arguments for
 *    use in the splineGrad functions.
 *
 */

short phaseGrad2est(double *ppv,double *u2,double *v2,int numRows,int numCols,
                    short horizSym,short vertSym) {
  double        *cosPhi;                /* cosine of the phase */
  double        *sinPhi;                /* sine of the phase */
  double        *u2cosPhi;              /* 2nd horizontal deriv. of cos(phi) */
  double        *v2cosPhi;              /* 2nd vertical deriv. of cos(phi) */
  int            pixCount;              /* pixel loop counter */
  int            numPix;                /* number of pixels in the image */

  numPix = numRows*numCols;

  /* allocate working arrays */
  if((cosPhi = (double *)malloc(numPix*sizeof(double))) == NULL) {
    printf("\nphaseGrad2est: Free store exhausted.\n");
    return -1;
  }
  if((sinPhi = (double *)malloc(numPix*sizeof(double))) == NULL) {
    printf("\nphaseGrad2est: Free store exhausted.\n");
    free(cosPhi);
    return -1;
  }
  if((u2cosPhi = (double *)malloc(numPix*sizeof(double))) == NULL) {
    printf("\nphaseGrad2est: Free store exhausted.\n");
    free(cosPhi);
    free(sinPhi);
    return -1;
  }
  if((v2cosPhi = (double *)malloc(numPix*sizeof(double))) == NULL) {
    printf("\nphaseGrad2est: Free store exhausted.\n");
    free(cosPhi);
    free(sinPhi);
    free(u2cosPhi);
    return -1;
  }

  /* calculate the sine and cosine of the phase */
  for(pixCount = 0; pixCount < numPix; pixCount++) {
    cosPhi[pixCount] = cos(ppv[pixCount]);
    sinPhi[pixCount] = sin(ppv[pixCount]);
  }

  /* compute the second phase derivatives */
  if(splineGrad2(sinPhi,u2,v2,numRows,numCols,horizSym,vertSym) != 0) {
    printf("\nphaseGrad2est: Sine second derivative failed.\n");
    free(cosPhi);
    free(sinPhi);
    free(u2cosPhi);
    free(v2cosPhi);
    return -1;
  }
  if(splineGrad2(cosPhi,u2cosPhi,v2cosPhi,numRows,numCols,horizSym,vertSym)
     != 0) {
    printf("\nphaseGrad2est: Cosine second derivative failed.\n");
    free(cosPhi);
    free(sinPhi);
    free(u2cosPhi);
    free(v2cosPhi);
    return -1;
  }
  for(pixCount = 0; pixCount < numPix; pixCount++) {
    u2[pixCount] = u2[pixCount]*cosPhi[pixCount] -
                   sinPhi[pixCount]*u2cosPhi[pixCount];
    v2[pixCount] = v2[pixCount]*cosPhi[pixCount] -
                   sinPhi[pixCount]*v2cosPhi[pixCount];
  }

  free(cosPhi);
  free(sinPhi);
  free(u2cosPhi);
  free(v2cosPhi);
  return 0;
}
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: splineGrad2
 *
 *  Function computes the horizontal and vertical second derivatives of the
 *    cubic spline whose knot values are given by img. numRows and numCols are
 *    the number of image rows and columns, respectively, and the second
 *    derivative images are returned in dx2 and dy2. Returns -1 on error.
 *
 *  4/25/2005  ras
 *
 *  REVISED 8/12/2005 ras: added horizontal and vertical mirror antisymmetry and
 *    no symmetry boundary conditions for the input data.
 *
 *    Symmetry		horizSym/vertSym
 *    --------		----------------
 *    None		0
 *    Mirror		positive
 *    Mirror Anti	negative
 *
 */

short splineGrad2(double *img,double *dx2,double *dy2,int numRows,int numCols,
                  short horizSym,short vertSym) {
  double       *colData;            /* temp column data array */
  double       *derivData;          /* temp derivative data array */
  int           row;                /* row counter */
  int           col;                /* col counter */
  
  /* check input size */
  if(numRows < 2) {
    printf("\nsplineGrad2: Insufficient number of rows.\n");
    return -1;
  }
  if(numCols < 2) {
    printf("\nsplineGrad2: Insufficient number of columns.\n");
    return -1;
  }
  
  /* allocate temporary column and derivative arrays */
  if((colData = (double *)malloc(numRows*sizeof(double))) == NULL) {
    printf("\nsplineGrad2: Free store exhausted.\n");
    return -1;
  }
  if((derivData = (double *)malloc(numRows*sizeof(double))) == NULL) {
    printf("\nsplineGrad2: Free store exhausted.\n");
    free(colData);
    return -1;
  }
  
  /* compute the second horizontal derivatives */
  for(row = 0; row < numRows; row++) {
    if(splineDeriv2(img+row*numCols,dx2+row*numCols,numCols,horizSym) != 0) {
      printf("\nsplineGrad2: Second horizontal derivative failed on row %d.\n",
             row);
      free(colData);
      free(derivData);
      return -1;
    }
  }
  
  /* compute the second vertical derivatives */
  for(col = 0; col < numCols; col++) {
    for(row = 0; row < numRows; row++) {
      colData[row] = img[row*numCols + col];
    }
    if(splineDeriv2(colData,derivData,numRows,vertSym) != 0) {
      printf("\nsplineGrad2: Second vertical derivative failed on column %d.\n",
             col);
      free(colData);
      free(derivData);
      return -1;
    }
    for(row = 0; row < numRows; row++) {
      dy2[row*numCols + col] = derivData[row];
    }
  }
  
  free(colData);
  free(derivData);
  return 0;
}
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: splineDeriv2
 *
 *  Function computes the second derivative of the cubic spline whose knot 
 *    values are given in array x of length len. The second derivative values at
 *    the knots are returned in xdotdot. See
 *      (1) "B-Spline Signal Processing: Part I..." pp. 824,826-827
 *      (2) "B-Spline Signal Processing: Part II..." pp. 835
 *    Returns -1 on error.
 *
 *  4/25/2005  ras
 *
 *  REVISED 8/12/2005 ras: added mirror antisymmetry and no symmetry boundary
 *    conditions for the input data. See pages 117 through 123 of thesis
 *    notebook.
 *
 *    Symmetry		sym
 *    --------		---
 *    None		0
 *    Mirror		positive
 *    Mirror Anti	negative
 *
 */

short splineDeriv2(double *x,double *xdotdot,int len,short sym) {
  double       *c3;                 /* cubic B-spline coefficients */
  int           i;                  /* loop counter */
  
  /* check input data length */
  if(len < 2) {
    printf("\nsplineDeriv2: Insufficient data length.\n");
    return -1;
  }
  
  /* allocate cubic B-spline coefficient array */
  if((c3 = (double *)malloc(len*sizeof(double))) == NULL) {
    printf("\nsplineDeriv2: Free store exhausted.\n");
    return -1;
  }
  
  /* compute cubic B-spline coefficients */
  if(symExpFilt(x,c3,len,sym,-2 + sqrt(3),6) != 0) {
    printf("\nsplineDeriv2: Computation of spline coefficients failed.\n");
    free(c3);
    return -1;
  }
  
  /* compute the second derivative samples */
  if(sym == 0) {                    /* no symmetry */
    xdotdot[0] = c3[1] - 2*c3[0] + c3[len-1];
    xdotdot[len-1] = c3[0] - 2*c3[len-1] + c3[len-2];
  } else {
    if(sym > 0) {                   /* mirror symmetry */
      xdotdot[0] = 2*(c3[1] - c3[0]);
      xdotdot[len-1] = 2*(c3[len-2] - c3[len-1]);
    } else {                        /* mirror antisymmetry */
      xdotdot[0] = 0;
      xdotdot[len-1] = 0;
    }
  }
  for(i = 1; i < len-1; i++) {
    xdotdot[i] = c3[i+1] - 2*c3[i] + c3[i-1];
  }
  
  free(c3);
  return 0;
}
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: splineGrad
 *
 *  Function computes the horizontal and vertical gradients of the cubic
 *    spline whose knot values are given by img. numRows and numCols are
 *    the number of image rows and columns, respectively, and the
 *    gradient images are returned in dx and dy. Returns -1 on error.
 *
 *  5/28/2004  ras
 *
 *  REVISED 4/25/2005 ras: added error messages for derivative computation
 *    failure.
 *
 *  REVISED 8/12/2005 ras: added horizontal and vertical mirror antisymmetry and
 *    no symmetry boundary conditions for the input data.
 *
 *    Symmetry		horizSym/vertSym
 *    --------		----------------
 *    None		0
 *    Mirror		positive
 *    Mirror Anti	negative
 *
 */

short splineGrad(double *img,double *dx,double *dy,int numRows,int numCols,
                 short horizSym,short vertSym) {
  double       *colData;            /* temp column data array */
  double       *gradData;           /* temp gradient data array */
  int           row;                /* row counter */
  int           col;                /* col counter */
  
  /* check input size */
  if(numRows < 2) {
    printf("\nsplineGrad: Insufficient number of rows.\n");
    return -1;
  }
  if(numCols < 2) {
    printf("\nsplineGrad: Insufficient number of columns.\n");
    return -1;
  }
  
  /* allocate temporary column and gradient arrays */
  if((colData = (double *)malloc(numRows*sizeof(double))) == NULL) {
    printf("\nsplineGrad: Free store exhausted.\n");
    return -1;
  }
  if((gradData = (double *)malloc(numRows*sizeof(double))) == NULL) {
    printf("\nsplineGrad: Free store exhausted.\n");
    free(colData);
    return -1;
  }
  
  /* compute the horizontal gradients */
  for(row = 0; row < numRows; row++) {
    if(splineDeriv(img+row*numCols,dx+row*numCols,numCols,horizSym) != 0) {
      printf("\nsplineGrad: Horizontal gradient failed on row %d.\n",row);
      free(colData);
      free(gradData);
      return -1;
    }
  }
  
  /* compute the vertical gradients */
  for(col = 0; col < numCols; col++) {
    for(row = 0; row < numRows; row++) {
      colData[row] = img[row*numCols + col];
    }
    if(splineDeriv(colData,gradData,numRows,vertSym) != 0) {
      printf("\nsplineGrad: Vertical gradient failed on column %d.\n",col);
      free(colData);
      free(gradData);
      return -1;
    }
    for(row = 0; row < numRows; row++) {
      dy[row*numCols + col] = gradData[row];
    }
  }
  
  free(colData);
  free(gradData);
  return 0;
}
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: splineDeriv
 *
 *  Function computes the derivative of the cubic spline whose knot 
 *    values are given in array x of length len. The derivative values at the
 *    knots are returned in xdot. See
 *      (1) "B-Spline Signal Processing: Part I..." pp. 824,826-827
 *      (2) "B-Spline Signal Processing: Part II..." pp. 835
 *    Returns -1 on error.
 *
 *  5/28/2004  ras
 *
 *  REVISED 4/25/2005 ras: added error message for failure to compute spline
 *    coefficients.
 *
 *  REVISED 8/12/2005 ras: added mirror antisymmetry and no symmetry boundary
 *    conditions for the input data. See pages 117 through 123 of thesis
 *    notebook.
 *
 *    Symmetry		sym
 *    --------		---
 *    None		0
 *    Mirror		positive
 *    Mirror Anti	negative
 *
 */

short splineDeriv(double *x,double *xdot,int len,short sym) {
  double       *c3;                 /* cubic B-spline coefficients */
  int           i;                  /* loop counter */
  
  /* check input data length */
  if(len < 2) {
    printf("\nsplineDeriv: Insufficient data length.\n");
    return -1;
  }
  
  /* allocate cubic B-spline coefficient array */
  if((c3 = (double *)malloc(len*sizeof(double))) == NULL) {
    printf("\nsplineDeriv: Free store exhausted.\n");
    return -1;
  }
  
  /* compute cubic B-spline coefficients */
  if(symExpFilt(x,c3,len,sym,-2 + sqrt(3),6) != 0) {
    printf("\nsplineDeriv: Computation of spline coefficients failed.\n");
    free(c3);
    return -1;
  }
  
  /* compute the derivative samples */
  if(sym == 0) {                    /* no symmetry */
    xdot[0] = (c3[1] - c3[len-1])/2;
    xdot[len-1] = (c3[0] - c3[len-2])/2;
  } else {
    if(sym > 0) {                   /* mirror symmetry */
      xdot[0] = 0;
      xdot[len-1] = 0;
    } else {                        /* mirror antisymmetry */
      xdot[0] = c3[1];
      xdot[len-1] = -c3[len-2];
    }
  }
  for(i = 1; i < len-1; i++) {
    xdot[i] = (c3[i+1] - c3[i-1])/2;
  }
  
  free(c3);
  return 0;
}
/* ----------------------------------------------------------------------------*/





/* ----------------------------------------------------------------------------*/
/*
 *  Function: symExpFilt
 *
 *  Function filters the input data x of length len with the second-order
 *    symmetrical exponential filter described by pole and gain parameters.
 *    This implementation is described in
 *      (1) "Splines - A Perfect Fit for Signal and Image Processing" pp. 26
 *      (2) "B-Spline Signal Processing: Part II..." pp. 835-836
 *    The filter output is stored at pointer cMinus. Returns -1 on error.
 *
 *  5/28/2004  ras
 *
 *  REVISED 4/4/2005 ras: simplified index calculation in anticausal filter
 *    loop.
 *
 *  REVISED 8/12/2005 ras: added mirror antisymmetry and no symmetry boundary
 *    conditions for the input data. See pages 117 through 123 of thesis
 *    notebook.
 *
 *    Symmetry		sym
 *    --------		---
 *    None		0
 *    Mirror		positive
 *    Mirror Anti	negative
 *
 */

short symExpFilt(double *x,double *cMinus,int len,short sym,double pole,
                 double gain) {
  double       *cPlus;              /* causal filter output */
  double        polePow;            /* aggregate powers of pole */
  int           i;                  /* loop counter */
  
  
  /* check input data length */
  if(len < 2) {
    printf("\nsymExpFilt: Insufficient data length.\n");
    return -1;
  }
  
  /* check that mirror antisymmetric (odd) input has zero-valued end-points */
  if((sym < 0) &&
     ((fabs(x[0]) >= DBL_EPSILON) || (fabs(x[len-1]) >= DBL_EPSILON))) {
    printf("\nsymExpFilt: Incompatible symmetry for data with non-zero-valued "
           "end-points.\n");
    return -1;
  }

  /* allocate causal filter output array */
  if((cPlus = (double *)malloc(len*sizeof(double))) == NULL) {
    printf("\nsymExpFilt: Free store exhausted.\n");
    return -1;
  }
  
  /* initialize causal filter */
  cPlus[0] = x[0];
  polePow = pole;
  if(sym == 0) {                    /* no symmetry */
    for(i = len-1; i > 0; i--) {
      cPlus[0] += x[i]*polePow;
      polePow *= pole;
    }
    cPlus[0] /= 1 - polePow;
  } else {                          /* mirror symmetry/antisymmetry */
    for(i = 1; i < len; i++) {
      cPlus[0] += x[i]*polePow;
      polePow *= pole;
    }
    if(sym > 0) {                   /* mirror symmetry */
      for(i = len-2; i > 0; i--) {
        cPlus[0] += x[i]*polePow;
        polePow *= pole;
      }
      cPlus[0] /= 1 - polePow;
    } else {                        /* mirror antisymmetry */
      for(i = len-2; i > 0; i--) {
        cPlus[0] -= x[i]*polePow;
        polePow *= pole;
      }
      cPlus[0] /= polePow - 1;
    }
  }

  /* apply causal recursion */
  for(i = 1; i < len; i++) {
    cPlus[i] = x[i] + pole*cPlus[i-1];
  }
  
  /* initialize anticausal filter */
  if(sym == 0) {                    /* no symmetry */
    cMinus[len-1] = x[0];
    polePow = pole;
    for(i = 1; i < len; i++) {
      cMinus[len-1] += x[i]*polePow;
      polePow *= pole;
    }
    cMinus[len-1] /= 1 - polePow;
    cMinus[len-1] = cPlus[len-1] + pole*cMinus[len-1];
  } else {
    if(sym > 0) {                   /* mirror symmetry */
      cMinus[len-1] = cPlus[len-1] + pole*cPlus[len-2];
    } else {                        /* mirror antisymmetry */
      cMinus[len-1] = cPlus[len-1] - pole*cPlus[len-2];
    }
  }
  cMinus[len-1] *= -gain*pole/(1 - pole*pole);

  /* apply anticausal recursion */
  for(i = len-2; i >= 0; i--) {
    cMinus[i] = pole*(cMinus[i+1] - gain*cPlus[i]);
  }
  
  /* equivalent, less efficient index calculation */
  /*
  for(i = 2; i <= len; i++) {
    cMinus[len-i] = pole*(cMinus[len+1-i] - gain*cPlus[len-i]);
  }
  */
  
  free(cPlus);
  return 0;
}
/* ----------------------------------------------------------------------------*/





