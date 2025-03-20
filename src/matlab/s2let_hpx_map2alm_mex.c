// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/**
 * MATLAB interface: s2let_hpx_map2alm_mex.
 * This function for internal use only.
 *
 * Usage:
 *   alm = s2let_hpx_map2alm_mex(f, nside, L);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int i, nside, L, f_m, f_n;
  complex double *flm = NULL;
  double *f_real = NULL, *flm_real = NULL, *flm_imag = NULL, *f_r = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=3) {
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_mex:InvalidInput:nrhs",
          "Require six inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }
  // Parse input dataset f
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_real = mxGetPr(prhs[iin]);
  f_r = (double*)malloc(f_m * f_n * sizeof(double));
  for (i=0; i<f_m*f_n; i++)
    f_r[i] = f_real[i];

  // Parse HEALPIX parameter nside
  iin = 1;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_mex:InvalidInput:healpixParameter",
          "HEALPIX parameter nside must be integer.");
  }
  nside = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)nside || nside <= 1)
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_mex:InvalidInput:healpixParameter",
          "Healpix parameter nside must be positive integer greater than 2");

  if( f_m*f_n != 12*nside*nside )
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_mex:InvalidInput:LbandLimit",
          "nside must correspond to the sampling scheme, i.e. f = 12*nside*nside samples.");

  // Parse harmonic band-limit L
  iin = 2;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  // Perform harmonic transform
  s2let_allocate_lm(&flm, L);
  s2let_hpx_map2alm_real(flm, f_r, nside, L);

  // Output flm's
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(L * L, 1, mxCOMPLEX);
  flm_real = mxGetPr(plhs[iout]);
  flm_imag = mxGetPi(plhs[iout]);
  for (i=0; i<L*L; i++){
    flm_real[i] = creal( flm[i] );
    flm_imag[i] = cimag( flm[i] );
  }

  free(flm);
  free(f_r);


}
