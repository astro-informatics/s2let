// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <ssht/ssht.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "s2let/s2let.h"

void s2let_lm_random_flm_real_sigma(
    complex double *flm, int L, int seed, double sigmanoise) {
  int el, m, msign, i, i_op;
  for (el = 0; el < L; el++) {
    m = 0;
    i = el * el + el + m;
    flm[i] = (2.0 * ran2_dp(seed) - 1.0);
    for (m = 1; m <= el; m++) {
      i = el * el + el + m;
      flm[i] = sigmanoise * (2.0 * ran2_dp(seed) - 1.0) +
               I * sigmanoise * (2.0 * ran2_dp(seed) - 1.0);
      i_op = el * el + el - m;
      msign = m & 1;
      msign = 1 - msign - msign; // (-1)^m
      flm[i_op] = msign * conj(flm[i]);
    }
  }
}

double needletpower(double *wav_lm, int L) {
  int i;
  double totalpower = 0;
  for (i = 0; i < L; i++)
    totalpower += pow(wav_lm[i], 2.0);
  return totalpower;
}

/*!
 * PROGRAM : s2let_axisym_denoising_demo
 * COMMAND : bin/s2let_axisym_denoising_demo
 * ARGUMENTS : none
 */
int main(int argc, char *argv[]) {
  s2let_parameters_t parameters = {};

  const int seed = (int)(10000.0 * (double)clock() / (double)CLOCKS_PER_SEC);
  // PARAMETERS
  const double SNR_in = 10.0; // Input SNR
  const int nsigma = 3;       // Number of sigmas for hard thresholding
  const int multires = 1;     // Multiresolution flag
  const double B = 2;         // Wavelet parameters
  const int J_min = 0;        // First wavelet scale to use

  char outfile[100];
  double *f, *noise, *g, *g_wav, *g_scal, *wav_lm, *scal_lm, *f_denois,
      *remaining_noise;
  complex double *noise_lm;

  parameters.B = B;
  parameters.J_min = J_min;

  printf("--------------------------------------------------\n");
  printf(" S2LET library : denoising example\n");
  printf(" Earth tomography signal, MW sampling\n");
  printf("--------------------------------------------------\n");

  char file[100] = "data/earth_tomo_mw_128.fits";
  printf(" Reading file %s\n", file);
  const int L = s2let_fits_mw_read_bandlimit(file);
  parameters.L = L;
  printf(" - Detected bandlimit L = %i\n", L);
  int J = s2let_j_max(&parameters);
  printf(" Parameters for wavelet denoising :\n");
  s2let_switch_wavtype(1);
  printf(" - Input SNR : %f\n", SNR_in);
  printf(" - Sigma threshold : %i\n", nsigma);
  printf(" - Multiresolution flag : %i\n", multires);
  printf(" - Wavelet parameter : %f\n", B);
  printf(" - Total number of wavelets : %i\n", J);
  printf(" - First wavelet scale to be used : %i\n", J_min);

  s2let_allocate_mw_real(&f, L);
  s2let_fits_mw_read_map(f, file, L); // Read MW map from file
  printf(" File successfully read from file\n");

  // Compute noise standard deviation and generate noise
  double sigmanoise = sqrt(pow(10.0, -SNR_in / 10.0) * s2let_mw_power_real(f, L));
  printf(" - Std dev of the noise (to match SNR) = %f\n", sigmanoise);
  s2let_allocate_lm(&noise_lm, L);
  s2let_lm_random_flm_real_sigma(noise_lm, L, seed, sigmanoise);
  double SNR_actual =
      10.0 * log10(s2let_mw_power_real(f, L) / s2let_lm_power(noise_lm, L));
  printf(" - Actual (realised) SNR = %f\n", SNR_actual);

  // Add noise to the signal in real space
  printf(" Contaminating the signal with this noise...");
  fflush(NULL);
  s2let_allocate_mw_real(&noise, L);
  s2let_mw_alm2map_real(noise, noise_lm, L);
  s2let_allocate_mw_real(&g, L);
  int i, j;
  for (i = 0; i < L * (2 * L - 1); i++)
    g[i] = f[i] + noise[i];
  printf(" done\n");

  printf(" Performing wavelet decomposition...");
  fflush(NULL);
  // Perform wavelet analysis from scratch with all signals given as MW maps
  if (multires) {
    s2let_transform_axisym_allocate_mw_f_wav_multires_real(
        &g_wav, &g_scal, &parameters);
    s2let_transform_axisym_wav_analysis_mw_multires_real(g_wav, g_scal, g, &parameters);
  } else {
    s2let_transform_axisym_allocate_mw_f_wav_real(&g_wav, &g_scal, &parameters);
    s2let_transform_axisym_wav_analysis_mw_real(g_wav, g_scal, g, &parameters);
  }
  printf(" done\n");

  // Compute simple threshold for needlet coefficients based on noise model
  printf(" Construct the threshold rule for the Gaussian noise\n");
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);
  double *treshold = (double *)calloc((J - J_min + 1), sizeof(double));
  for (j = J_min; j <= J; j++)
    treshold[j - J_min] = sigmanoise * nsigma * sqrt(needletpower(wav_lm + j * L, L));

  printf(" Hard thresholding the wavelets...");
  fflush(NULL);
  s2let_allocate_mw_real(&f_denois, L);
  if (multires) {
    s2let_transform_axisym_wav_hardthreshold_multires_real(
        g_wav, treshold, &parameters);
    s2let_transform_axisym_wav_synthesis_mw_multires_real(
        f_denois, g_wav, g_scal, &parameters);
  } else {
    s2let_transform_axisym_wav_hardthreshold_real(g_wav, treshold, &parameters);
    s2let_transform_axisym_wav_synthesis_mw_real(f_denois, g_wav, g_scal, &parameters);
  }
  printf(" done\n");

  // Remaining noise
  s2let_allocate_mw_real(&remaining_noise, L);
  for (i = 0; i < L * (2 * L - 1); i++)
    remaining_noise[i] = f_denois[i] - f[i];

  // SNR after denoising
  double SNR_denoised =
      10.0 * log10(s2let_mw_power_real(f, L) / s2let_mw_power_real(remaining_noise, L));
  printf(" -> SNR before denoising = %f\n", SNR_actual);
  printf(" -> SNR after denoising  = %f\n", SNR_denoised);

  // Finally write the denoised signal
  printf(" Write output files\n");
  sprintf(outfile, "%s%s%s", "data/earth_tomo_mw_128", "_noisy", ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile);                        // In case the file exists
  s2let_fits_mw_write_map(outfile, g, L); // Now write the map to fits file
  char params[100];
  sprintf(params, "%d%s%d%s%d", L, "_", B, "_", J_min);
  sprintf(outfile, "%s%s%s", "data/earth_tomo_mw_128", "_denoised", ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile);                               // In case the file exists
  s2let_fits_mw_write_map(outfile, f_denois, L); // Now write the map to fits file

  free(f);
  free(noise);
  free(g);
  free(g_wav);
  free(g_scal);
  free(scal_lm);
  free(f_denois);
  free(remaining_noise);
  free(noise_lm);

  printf("--------------------------------------------------\n");
  return 0;
}
