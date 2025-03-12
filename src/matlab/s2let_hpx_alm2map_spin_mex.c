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
 *   f = s2let_hpx_alm2map_spin_mex(almQ, almU, nside, L, spin);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int i, nside, spin, L, f_m, f_n;
  complex double *flmQ = NULL, *flmU = NULL;
  double *flm_imag = NULL, *flm_real = NULL, *f_real = NULL, *fQ_r = NULL, *fU_r = NULL;
  int iin = 0, iout = 0;


  // Parse input 
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  flm_real = mxGetPr(prhs[iin]);
  flm_imag = mxGetPi(prhs[iin]);
  flmQ = (complex double*)malloc( f_m*f_n * sizeof(complex double));
    for(i=0; i<f_n*f_m; i++)
      flmQ[ i ] = flm_real[ i ] + I * flm_imag[ i ] ;

  // Parse input 
  iin = 1;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  flm_real = mxGetPr(prhs[iin]);
  flm_imag = mxGetPi(prhs[iin]);
  flmU = (complex double*)malloc( f_m*f_n * sizeof(complex double));
    for(i=0; i<f_n*f_m; i++)
      flmU[ i ] = flm_real[ i ] + I * flm_imag[ i ] ;

  // Parse HEALPIX parameter nside
  iin = 2;
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
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_hpx_alm2map_spin_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_hpx_alm2map_spin_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  // Parse spin
  iin = 4;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_hpx_alm2map_spin_mex:InvalidInput:LbandLimit",
          "spin must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);

  // Perform harmonic transform 
  s2let_hpx_allocate_real(&fQ_r, nside);
  s2let_hpx_allocate_real(&fU_r, nside);
  s2let_hpx_alm2map_spin_real(fQ_r, fU_r, flmQ, flmU, nside, L, spin);
   
  // Output function f
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, 12 * nside * nside, mxREAL);
  f_real = mxGetPr(plhs[iout]);
  for (i=0; i < 12 * nside * nside; i++)
    f_real[i] = creal(fQ_r[i]);

  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(1, 12 * nside * nside, mxREAL);
  f_real = mxGetPr(plhs[iout]);
  for (i=0; i < 12 * nside * nside; i++)
    f_real[i] = creal(fU_r[i]);

  free(fQ_r);
  free(flmQ);
  free(fU_r);
  free(flmU);

}
