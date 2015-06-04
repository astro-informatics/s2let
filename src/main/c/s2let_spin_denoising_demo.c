// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void s2let_lm_random_flm_sigma(complex double *flm, int L, int seed, double sigmanoise) {
  int i;
  for (i = 0; i < L*L; ++i)
    flm[i] = sigmanoise * (2.0*ran2_dp(seed) - 1.0) + I * sigmanoise * (2.0*ran2_dp(seed) - 1.0);
}

double waveletpower(complex double *wav_lm, int L){
  int i;
  double totalpower = 0;
  for(i = 0; i < L*L; i++)
    totalpower += wav_lm[i] * conj(wav_lm[i]);
  return totalpower;
}

void hard_threshold(
    complex double *g_wav,
    const double *threshold,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    int J = s2let_j_max(parameters);
    int i, j, offset = 0;
    for(j = J_min; j <= J; j++){
        int bl = parameters->upsample ? L : MIN(s2let_bandlimit(j, parameters), L);
        for(i = 0; i < N*bl*(2*bl-1); i++){
            if( cabs(g_wav[offset + i]) < threshold[j-J_min] )
                g_wav[offset + i] = 0;
        }
        offset += N*bl*(2*bl-1);
    }
}

/*!
 * PROGRAM : s2let_denoising_demo
 * COMMAND : bin/s2let_denoising_demo
 * ARGUMENTS : none
 */
int main(int argc, char *argv[])
{
  int i, j;
  s2let_parameters_t parameters = {};

  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
  // PARAMETERS
  const double SNR_in = 10.0;  // Input SNR
  const int nsigma = 3;   // Number of sigmas for hard thresholding
  const int multires = 1; // Multiresolution flag
  const int B = 2;        // Wavelet parameters
  const int N = 4;        // Azimuthal band-limit
  const int J_min = 0;    // First wavelet scale to use
  const int spin = 2;     // Spin number

  char outfile[100];
  complex double *flm, *f, *noise, *g, *g_wav, *g_scal, *f_denoised, *remaining_noise;
  double *f_hpx_r, *f_hpx_i,
         *f_r, *f_i,
         *noise_r, *noise_i,
         *f_denoised_r, *f_denoised_i,
         *g_r, *g_i,
         *scal_l;
  complex double *noise_lm, *wav_lm;

  parameters.B = B;
  parameters.J_min = J_min;
  parameters.N = N;
  parameters.spin = spin;
  parameters.upsample = !multires;

  printf("--------------------------------------------------\n");
  printf(" S2LET library : denoising example\n");
  printf(" Earth tomography signal, MW sampling\n");
  printf("--------------------------------------------------\n");

  // Commented out code uses a real map and converts it to
  // a spin map

  //char file[100] = "data/earth_tomo_mw_128_rot.fits";
  char fileQ[100] = "data/wmap_mcmc_base_k_synch_stk_q_9yr_v5_rot_denoised.fits";
  char fileU[100] = "data/wmap_mcmc_base_k_synch_stk_u_9yr_v5_rot_denoised.fits";
  //printf(" Reading file %s\n", file);
  printf(" Reading file %s\n", fileQ);
  printf(" Reading file %s\n", fileU);
  const int L = s2let_fits_mw_read_bandlimit(fileQ);
  //const int nside = s2let_fits_hpx_read_nside(fileQ);
  //int L = MIN(192,3*nside);
  parameters.L = L;
  printf(" - Detected bandlimit L = %i\n",L);
  int J = s2let_j_max(&parameters);
  printf(" Parameters for wavelet denoising :\n");
  s2let_switch_wavtype(1);
  printf(" - Input SNR : %f\n",SNR_in);
  printf(" - Sigma threshold : %i\n", nsigma);
  printf(" - Multiresolution flag : %i\n", multires);
  printf(" - Wavelet parameter : %i\n", B);
  printf(" - Total number of wavelets : %i\n", J);
  printf(" - First wavelet scale to be used : %i\n", J_min);

  //s2let_allocate_mw_real(&f_r, L);
  //s2let_fits_mw_read_map(f_r, file, L); // Read MW map from file
  //s2let_hpx_allocate_real(&f_hpx_r, nside);
  //s2let_hpx_allocate_real(&f_hpx_i, nside);
  //s2let_hpx_read_map(f_hpx_r, fileQ, nside);
  //s2let_hpx_read_map(f_hpx_i, fileU, nside);
  //printf(" Map successfully read from file\n");

  // Force into spin signal
  // s2let_allocate_lm(&flm, L);
  // s2let_mw_map2alm_real(flm, f_r, L);
  // for (i = 0; i < spin*spin; ++i)
  //   flm[i] = 0;
  // s2let_allocate_mw(&f, L);
  // int verbosity = 0;
  // ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  // ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

  // Convert to MW sampling
  s2let_allocate_mw_real(&f_r, L);
  s2let_allocate_mw_real(&f_i, L);
  s2let_fits_mw_read_map(f_r, fileQ, L);
  s2let_fits_mw_read_map(f_i, fileU, L);
  //s2let_allocate_lm(&flm, L);
  //s2let_hpx_map2alm_real(flm, f_hpx_r, nside, L);
  //s2let_mw_alm2map_real(f_r, flm, L);
  //s2let_hpx_map2alm_real(flm, f_hpx_i, nside, L);
  //s2let_mw_alm2map_real(f_i, flm, L);

  s2let_allocate_mw(&f, L);
  for (i = 0; i < L*(2*L-1); ++i)
    f[i] = f_r[i] + I*f_i[i];
  //free(flm);
  //free(f_hpx_r);
  //free(f_hpx_i);

  // Compute noise standard deviation and generate noise
  double sigmanoise = sqrt(pow(10.0, -SNR_in/10.0) * s2let_mw_power(f, L));
  printf(" - Std dev of the noise (to match SNR) = %f\n", sigmanoise);
  s2let_allocate_lm(&noise_lm, L);
  s2let_lm_random_flm_sigma(noise_lm, L, seed, sigmanoise);
  double SNR_actual = 10.0 * log10( s2let_mw_power(f, L) / s2let_lm_power(noise_lm, L));
  printf(" - Actual (realised) SNR = %f\n", SNR_actual);

  // Add noise to the signal in real space
  printf(" Contaminating the signal with this noise...");fflush(NULL);
  s2let_allocate_mw(&noise, L);
  s2let_mw_alm2map(noise, noise_lm, L, 0);
  s2let_allocate_mw(&g, L);
  for (i = 0; i < L*(2*L-1); i++)
    g[i] = f[i] + noise[i];
  printf(" done\n");

  printf(" Performing wavelet decomposition...");fflush(NULL);
  // Perform wavelet analysis from scratch with all signals given as MW maps
  s2let_allocate_f_wav(&g_wav, &g_scal, &parameters);
  s2let_analysis_px2wav(g_wav, g_scal, g, &parameters);
  printf(" done\n");

  // Compute simple threshold for needlet coefficients based on noise model
  printf(" Construct the threshold rule for the Gaussian noise\n");
  s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, &parameters);
  s2let_tiling_wavelet(wav_lm, scal_l, &parameters);
  double *threshold = (double*)calloc((J-J_min+1), sizeof(double));
  for(j = J_min; j <= J; j++)
    threshold[j-J_min] = sigmanoise * nsigma * sqrt(waveletpower(wav_lm + j * L*L, L));

  printf(" Hard thresholding the wavelets...");fflush(NULL);
  s2let_allocate_mw(&f_denoised, L);
  hard_threshold(g_wav, threshold, &parameters);
  s2let_synthesis_wav2px(f_denoised, g_wav, g_scal, &parameters);
  printf(" done\n");

  // Remaining noise
  s2let_allocate_mw(&remaining_noise, L);
  for(i = 0; i < L*(2*L-1); i++)
    remaining_noise[i] = f_denoised[i] - f[i];

  // SNR after denoising
  double SNR_denoised = 10.0 * log10( s2let_mw_power(f, L) / s2let_mw_power(remaining_noise, L));
  printf(" -> SNR before denoising = %f\n", SNR_actual);
  printf(" -> SNR after denoising  = %f\n", SNR_denoised);

  // Finally write the denoised signal
  s2let_allocate_mw_real(&noise_r, L);
  s2let_allocate_mw_real(&noise_i, L);
  s2let_allocate_mw_real(&g_r, L);
  s2let_allocate_mw_real(&g_i, L);
  s2let_allocate_mw_real(&f_denoised_r, L);
  s2let_allocate_mw_real(&f_denoised_i, L);
  for (i = 0; i < L*(2*L-1); ++i)
  {
    noise_r[i] = creal(noise[i]);
    noise_i[i] = cimag(noise[i]);
    g_r[i] = creal(g[i]);
    g_i[i] = cimag(g[i]);
    f_denoised_r[i] = creal(f_denoised[i]);
    f_denoised_i[i] = cimag(f_denoised[i]);
  }
  printf(" Write output files\n");
  sprintf(outfile, "%s%s%s", "data/spin_signal_real", "_input" , ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile); // In case the file exists
  s2let_fits_mw_write_map(outfile, f_r, L); // Now write the map to fits file
  sprintf(outfile, "%s%s%s", "data/spin_signal_imag", "_input" , ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile); // In case the file exists
  s2let_fits_mw_write_map(outfile, f_i, L); // Now write the map to fits file

  sprintf(outfile, "%s%s%s", "data/spin_signal_real", "_noise" , ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile); // In case the file exists
  s2let_fits_mw_write_map(outfile, noise_r, L); // Now write the map to fits file
  sprintf(outfile, "%s%s%s", "data/spin_signal_imag", "_noise" , ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile); // In case the file exists
  s2let_fits_mw_write_map(outfile, noise_i, L); // Now write the map to fits file

  sprintf(outfile, "%s%s%s", "data/spin_signal_real", "_input_noise" , ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile); // In case the file exists
  s2let_fits_mw_write_map(outfile, g_r, L); // Now write the map to fits file
  sprintf(outfile, "%s%s%s", "data/spin_signal_imag", "_input_noise" , ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile); // In case the file exists
  s2let_fits_mw_write_map(outfile, g_i, L); // Now write the map to fits file

  sprintf(outfile, "%s%s%s", "data/spin_signal_real", "_denoised", ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile); // In case the file exists
  s2let_fits_mw_write_map(outfile, f_denoised_r, L); // Now write the map to fits file
  sprintf(outfile, "%s%s%s", "data/spin_signal_imag", "_denoised", ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile); // In case the file exists
  s2let_fits_mw_write_map(outfile, f_denoised_i, L); // Now write the map to fits file

  free(f);
  free(f_r);
  free(f_i);
  free(f_denoised);
  free(f_denoised_r);
  free(f_denoised_i);
  free(noise);
  free(noise_r);
  free(noise_i);
  free(g);
  free(g_r);
  free(g_i);
  free(g_wav);
  free(g_scal);
  free(scal_l);
  free(remaining_noise);
  free(noise_lm);

  printf("--------------------------------------------------\n");
  return 0;
}
