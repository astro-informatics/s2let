// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include <s2let_mex.h>
#include <string.h>
#include "mex.h"

/**
 * MATLAB interface: s2let_transform_synthesis_lm2wav_mex.
 * This function for internal use only.
 * Compute spin directional wavelet transform (synthesis)
 * with input in harmonic space and output in pixel space.
 *
 * Usage:
 *   flm = ...
 *        s2let_transform_synthesis_lm2wav_mex(f_wav, f_scal, B, L, J_min, N, spin, reality, upsample,
 *                                             spin_lowered, original_spin, sampling_scheme);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int i, j, L, J_min, N, spin, f_m, f_n, reality, upsample, original_spin;
  double B;
  char sampling_str[S2LET_STRING_LEN];
  s2let_sampling_t sampling_scheme;
  s2let_parameters_t parameters = {};
  double *f_wav_real, *f_scal_real, *flm_real;
  double *f_wav_imag, *f_scal_imag, *flm_imag;
  complex double *f_wav = NULL, *f_scal = NULL, *flm = NULL;
  double *f_wav_r = NULL, *f_scal_r = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=11) {
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:nrhs",
          "Require eleven inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse reality flag
  iin = 7;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:reality",
          "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

  // Parse multiresolution flag
  iin = 8;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:upsample",
          "Multiresolution flag must be logical.");
  upsample = mxIsLogicalScalarTrue(prhs[iin]);

  /* Parse sampling scheme method. */
  iin = 10;
  if( !mxIsChar(prhs[iin]) ) {
      mexErrMsgIdAndTxt("s2let_transform_synthesis_mw_mex:InvalidInput:samplingSchemeChar",
                        "Sampling scheme must be string.");
  }
  int len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
  if (len >= S2LET_STRING_LEN)
      mexErrMsgIdAndTxt("s2let_transform_synthesis_mw_mex:InvalidInput:samplingSchemeTooLong",
                        "Sampling scheme exceeds string length.");
  mxGetString(prhs[iin], sampling_str, len);

  if (strcmp(sampling_str, S2LET_SAMPLING_MW_STR) == 0)
      sampling_scheme = S2LET_SAMPLING_MW;
  else if (strcmp(sampling_str, S2LET_SAMPLING_MW_SS_STR) == 0)
      sampling_scheme = S2LET_SAMPLING_MW_SS;
  else
      mexErrMsgIdAndTxt("s2let_transform_synthesis_mw_mex:InvalidInput:samplingScheme",
                        "Invalid sampling scheme.");

  // Parse original spin
  iin = 9;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:spinloweredfrom",
          "SpinLoweredFrom must be integer.");
  }
  original_spin = (int)mxGetScalar(prhs[iin]);

  // Parse input wavelets f_wav
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_wav_real = mxGetPr(prhs[iin]);
  if(reality){
    f_wav_r = (double*)malloc( f_m*f_n * sizeof(double));
    for(j=0; j<f_m*f_n; j++)
      f_wav_r[ j ] = f_wav_real[ j ];
  }else{
    f_wav_imag = mxGetPi(prhs[iin]);
    f_wav = (complex double*)malloc( f_m*f_n * sizeof(complex double));
    for(j=0; j<f_n*f_m; j++)
      f_wav[ j ] = f_wav_real[ j ]
          + I * f_wav_imag[ j ] ;
  }

  // Parse input scaling function f_scal
  iin = 1;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_scal_real = mxGetPr(prhs[iin]);
  if(reality){
    f_scal_r = (double*)malloc( f_m*f_n * sizeof(double));
    for (i=0; i<f_m*f_n; i++)
      f_scal_r[i] = f_scal_real[i];
  }else{
    f_scal_imag = mxGetPi(prhs[iin]);
    f_scal = (complex double*)malloc( f_m*f_n * sizeof(complex double));
    for (i=0; i<f_m*f_n; i++)
      f_scal[i] = f_scal_real[i] + I * f_scal_imag[i];
  }

  // Parse wavelet parameter B
  iin = 2;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be integer.");
  }
  B = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)B || B <= 1)
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 3;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  // Parse first scale J_min
  iin = 4;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:Jmin",
          "First scale J_min must be integer.");
  }
  J_min = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)J_min || J_min < 0)
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:Jmin",
          "First scale J_min must be positive integer.");

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  // Compute ultimate scale J_max
  int J = s2let_j_max(&parameters);

  if( J_min > J+1 ) {
    mexErrMsgIdAndTxt("s2let_transform_synthesis_lm2wav_mex:InvalidInput:Jmin",
          "First scale J_min must be larger than that!");
  }

   // Parse azimuthal/directional band-limit N
  iin = 5;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:NbandLimit",
          "Azimuthal/directional band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);

    // Parse spin
  iin = 6;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:spin",
          "spin must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);

  parameters.N = N;
  parameters.spin = spin;
  parameters.upsample = upsample;
  parameters.original_spin = original_spin;
  parameters.reality = reality;
  parameters.sampling_scheme = sampling_scheme;

  // Perform wavelet transform in harmonic space and then reconstruction.
  s2let_allocate_lm(&flm, L);
  if(reality){
      s2let_synthesis_wav2lm_real(flm, f_wav_r, f_scal_r, &parameters);
  }else{
      s2let_synthesis_wav2lm(flm, f_wav, f_scal, &parameters);
  }

  // Output flm
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(L*L, 1, mxCOMPLEX);

  flm_real = mxGetPr(plhs[iout]);
  flm_imag = mxGetPi(plhs[iout]);
  for (i = 0; i < L*L; ++i)
  {
    flm_real[i] = creal( flm[i] );
    flm_imag[i] = cimag( flm[i] );
  }


  free(flm);
  if(reality){
    free(f_wav_r);
    free(f_scal_r);
  }else{
    free(f_wav);
    free(f_scal);
  }

}
