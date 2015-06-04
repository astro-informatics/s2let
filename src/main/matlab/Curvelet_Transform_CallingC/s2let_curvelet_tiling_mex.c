// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen
// Modified by Jennifer: s2let_tiling_curvelet instead of s2let_tiling_wavelet
//                        (add function s2let_tiling_curvelet also in s2let_tiling.c)

#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: s2let_curvelet_tiling_mex.
 * This function for internal use only.
 * Compute tiling in harmonic space for curvelets
 *
 * Usage:
 *   [psi_lm phi_l] = s2let_curvelet_tiling_mex(B, L, N, spin, J_min, SpinLowered, original_spin);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int B, L, J_min, spin, N, normalization, original_spin;
  s2let_parameters_t parameters = {};
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=7) {
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
  B = (int)mxGetScalar(prhs[iin]);
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

  // Parse normalization flag
  iin = 5;
  if( !mxIsLogicalScalar(prhs[iin])) {
    mexErrMsgIdAndTxt("s2let_tiling_mex:InvalidInput:SpinLowered",
          "SpinLowered flag must be logical.");
  }
  if (mxIsLogicalScalarTrue(prhs[iin]))
    normalization = S2LET_WAV_NORM_SPIN_LOWERED;
  else
    normalization = S2LET_WAV_NORM_DEFAULT;

  // Parse original spin
  iin = 6;
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
  parameters.normalization = normalization;
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
  s2let_tiling_curvelet_allocate(&psi_lm, &phi_l, &parameters);


  // Run S2LET function
  // Curvelets coefficients stored in 2D->1D array psi_lm, indexing with "ind = el .* el + el + m + 1".
  // where size of psi_lm = ind_pm *J+l
  // Scaling coefficients stored in 1D array phi_l with size L.
  s2let_tiling_curvelet(psi_lm, phi_l, &parameters);
    

  // Output psi_lm and phi_l
  double *psi_lm_out_real, *psi_lm_out_imag, *phi_l_out;
  int l, j, el, m;

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(L*L, J+1, mxCOMPLEX);
  psi_lm_out_real = mxGetPr(plhs[iout]);
  psi_lm_out_imag = mxGetPi(plhs[iout]);
  for (j = J_min; j <= J; ++j)
  {
     // int ind = spin*spin;
      for (el = ABS(spin); el < L; ++el)
      {
          m=el;
          
          int ind_pm= el*el+el+m;    //  ind_pm= ssht_sampling_elm2ind(ind, el, m);
          psi_lm_out_real[j*L*L + ind_pm] =  creal( psi_lm[j*L*L + ind_pm] );
          psi_lm_out_imag[j*L*L + ind_pm] =  cimag( psi_lm[j*L*L + ind_pm] );
                                                   
                                                   
          int ind_nm=el*el+el-m;     //ind_nm = ssht_sampling_elm2ind(ind,el, -m);
          psi_lm_out_real[j*L*L + ind_nm] =  creal( psi_lm[j*L*L + ind_nm] );
          psi_lm_out_imag[j*L*L + ind_nm] =  cimag( psi_lm[j*L*L + ind_nm] );
          
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
