// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: s2let_wavelet_tiling_mex.
 * This function for internal use only.
 * Compute tiling in harmonic space for directional wavelets.
 *
 * Usage:
 *   [psi_lm phi_l] = s2let_wavelet_tiling_mex(B, L, N, spin, J_min, SpinLowered, original_spin);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int L, J_min, spin, N, original_spin;
  double B;
  s2let_parameters_t parameters = {};
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=6) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:nrhs",
		      "Require seven inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidOutput:nlhs",
		      "Require two outputs.");
  }

  // Parse wavelet parameter B
  iin = 0;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be integer.");
  }
  B = (double)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)B || B <= 1)
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 1;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit L must be positive integer.");

  if( B >= L ) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be smaller than L!");
  }

  // Parse directional band-limit N
  iin = 2;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:LbandLimit",
          "Directional band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);

  // Parse spin
  iin = 3;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:LbandLimit",
          "spin must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);

  // Parse first scale J_min
  iin = 4;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:Jmin",
          "First scale J_min must be integer.");
  }
  J_min = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)J_min || J_min < 0)
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:Jmin",
          "First scale J_min must be positive integer.");

  // Parse original spin
  iin = 5;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:SpinLoweredFrom",
          "SpinLoweredFrom must be integer.");
  }
  original_spin = (int)mxGetScalar(prhs[iin]);


  parameters.L = L;
  parameters.B = B;
  parameters.J_min = J_min;
  parameters.N = N;
  parameters.spin = spin;
  parameters.original_spin = original_spin;

  // Compute ultimate scale J_max
  int J = s2let_j_max(&parameters);

  if( J_min > J+1 ) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:Jmin",
          "First scale J_min must be larger than that!");
  }

  // Allocate arrays
  complex double *psi_lm;
  double *phi_l;
  s2let_tiling_wavelet_allocate(&psi_lm, &phi_l, &parameters);

  // Run S2LET function
  s2let_tiling_wavelet(psi_lm, phi_l, &parameters);

  // Output psi_lm and phi_l
  double *psi_lm_out_real, *psi_lm_out_imag, *phi_l_out;
  int l, j, el, m;

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(L*L, J+1, mxCOMPLEX);
  psi_lm_out_real = mxGetPr(plhs[iout]);
  psi_lm_out_imag = mxGetPi(plhs[iout]);
  for (j = J_min; j <= J; ++j)
  {
      int ind = spin*spin;
      for (el = ABS(spin); el < L; ++el)
      {
          for (m = -el; m <= el; ++m)
          {
              psi_lm_out_real[j*L*L + ind] = creal( psi_lm[j*L*L + ind] );
              psi_lm_out_imag[j*L*L + ind] = cimag( psi_lm[j*L*L + ind] );
              ++ind;
          }
      }
  }

  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(1, L, mxREAL);
  phi_l_out = mxGetPr(plhs[iout]);
  for(l=0; l<L; l++)
    phi_l_out[l] = phi_l[l];

  free(phi_l);
  free(psi_lm);

}
