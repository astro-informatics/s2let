// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: s2let_hpx_axisym_synthesis.
 * This function for internal use only.
 *
 * Usage: 
 *   f = s2let_hpx_alm2map_mex(alm, nside, L);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int i, nside, L, f_m, f_n;
  complex double *flm = NULL;
  double *flm_imag = NULL, *flm_real = NULL, *f_real = NULL, *f_r = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=3) {
    mexErrMsgIdAndTxt("s2let_hpx_alm2map_mex:InvalidInput:nrhs",
          "Require seven inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("s2let_hpx_alm2map_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse input wavelets f_wav
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  flm_real = mxGetPr(prhs[iin]);
  flm_imag = mxGetPi(prhs[iin]);
  flm = (complex double*)malloc( f_m*f_n * sizeof(complex double));
    for(i=0; i<f_n*f_m; i++)
      flm[ i ] = flm_real[ i ] + I * flm_imag[ i ] ;

  // Parse HEALPIX parameter nside
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:healpixParameter",
          "HEALPIX parameter nside must be integer.");
  }
  nside = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)nside || nside <= 1)
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:healpixParameter",
          "Healpix parameter nside must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_hpx_alm2map_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_hpx_alm2map_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  // Perform harmonic transform 
  s2let_hpx_allocate_real(&f_r, nside);
  s2let_hpx_alm2map_real(f_r, flm, nside, L);
   
  // Output function f
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, 12 * nside * nside, mxREAL);
  f_real = mxGetPr(plhs[iout]);
  for (i=0; i < 12 * nside * nside; i++)
    f_real[i] = creal(f_r[i]);

  free(f_r);
  free(flm);

}
