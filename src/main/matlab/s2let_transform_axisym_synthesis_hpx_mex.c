// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: s2let_transform_axisym_synthesis_hpx.
 * This function for internal use only.
 * Compute axisymmetric wavelet transform (synthesis_hpx)
 * with output in pixel space. Format : HEALPIX maps.
 *
 * Usage:
 *   f = ...
 *        s2let_transform_axisym_synthesis_hpx_mex(f_wav, f_scal, nside, B, L, J_min);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int i, j, nside, L, J_min, f_m, f_n;
  double B;
  s2let_parameters_t parameters = {};
  double *f_wav_real, *f_scal_real, *f_real;
  double *f_wav_r = NULL, *f_scal_r = NULL, *f_r = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=6) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidInput:nrhs",
          "Require seven inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }


  // Parse input wavelets f_wav
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_wav_real = mxGetPr(prhs[iin]);
  f_wav_r = (double*)malloc( f_m*f_n * sizeof(double));
  for(j=0; j<f_m*f_n; j++)
    f_wav_r[ j ] = f_wav_real[ j ];

  // Parse input scaling function f_scal
  iin = 1;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_scal_real = mxGetPr(prhs[iin]);
  f_scal_r = (double*)malloc( f_m*f_n * sizeof(double));
  for (i=0; i<f_m*f_n; i++)
    f_scal_r[i] = f_scal_real[i];

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

  // Parse wavelet parameter B
  iin = 3;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be integer.");
  }
  B = (double)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)B || B <= 1)
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 4;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  // Parse first scale J_min
  iin = 5;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidInput:Jmin",
          "First scale J_min must be integer.");
  }
  J_min = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)J_min || J_min < 0)
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidInput:Jmin",
          "First scale J_min must be positive integer.");

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  // Compute ultimate scale J_max
  int J = s2let_j_max(&parameters);

  if( J_min > J+1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_synthesis_hpx_mex:InvalidInput:Jmin",
          "First scale J_min must be larger than that!");
  }


  // Perform wavelet transform in harmonic space and then reconstruction.
  s2let_hpx_allocate_real(&f_r, nside);
  s2let_transform_axisym_wav_synthesis_hpx_real(f_r, f_wav_r, f_scal_r, nside, &parameters);


  // Output function f
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, 12 * nside * nside, mxREAL);
  f_real = mxGetPr(plhs[iout]);
  for (i=0; i < 12 * nside * nside; i++)
    f_real[i] = creal(f_r[i]);

  free(f_r);
  free(f_wav_r);
  free(f_scal_r);

}
