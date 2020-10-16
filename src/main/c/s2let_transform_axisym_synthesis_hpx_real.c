// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "s2let/s2let.h"

/*!
 * PROGRAM : s2let_transform_axisym_hpx_synthesis_real
 * COMMAND : bin/s2let_transform_axisym_hpx_synthesis_real file B J_min L
 * ARGUMENTS :
 * - file : input root for the healpix wavelet map
 * - B : wavelet parameter
 * - J_min : first wavelet scale to use
 * - L : bandlimit for the decomposition
 * OUTPUT : fits files containing the reconstructed healpix map
 */
int main(int argc, char *argv[]) {
  printf("--------------------------------------------------\n");
  printf("S2LET library : axisymmetric wavelet transform\n");
  printf("Real signal, HEALPIX sampling\n");
  printf("--------------------------------------------------\n");

  char fileroot[100];
  int L, B, J_min;
  if (sscanf(argv[1], "%s", fileroot) != 1)
    exit(-2);
  if (sscanf(argv[2], "%i", &B) != 1)
    exit(-2);
  if (sscanf(argv[3], "%i", &J_min) != 1)
    exit(-2);
  if (sscanf(argv[4], "%i", &L) != 1)
    exit(-2);

  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  printf("Parameters for wavelet transform :\n");
  printf("- Wavelet parameter : %i\n", L);
  int J = s2let_j_max(&parameters);
  printf("- Wavelet parameter : %i\n", B);
  printf("- Total number of wavelets : %i\n", J);
  printf("- First wavelet scale to be used : %i\n", J_min);

  char params[100];
  char file[100];
  sprintf(params, "%d%s%d%s%d", L, "_", B, "_", J_min);
  int j, offset = 0;
  printf("File root = %s\n", fileroot);

  // Init read procedure, read nside
  sprintf(file, "%s%s%s%s", fileroot, "_scal_", params, ".fits");
  printf("- Infile_scal = %s\n", file);
  const int nside = s2let_fits_hpx_read_nside(file);
  printf("- Detected bandlimit nside = %i\n", nside);
  // Allocate memory for wavelets
  double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_hpx_f_wav_real(&f_wav, &f_scal, nside, &parameters);
  // Read the scaling function
  s2let_hpx_read_map(f_scal, file, nside); // Now write the map to fits file
  // Read the wavelets
  for (j = J_min; j <= J; j++) {
    sprintf(file, "%s%s%s%s%d%s", fileroot, "_wav_", params, "_", j, ".fits");
    printf("- Infile_wav[j=%i] = %s\n", j, file);
    s2let_hpx_read_map(f_wav + offset, file, nside); // Now write the map to fits file
    offset += 12 * nside * nside;                    // Go to the next wavelet
  }

  // Allocate memory for reconstruction
  double *f = (double *)calloc(12 * nside * nside, sizeof(double));
  printf("File successfully read from file\n");

  printf("Performing wavelet reconstruction...");
  fflush(NULL);
  s2let_transform_axisym_wav_synthesis_hpx_real(f, f_wav, f_scal, nside, &parameters);
  printf("done\n");

  // Output the wavelets to FITS files
  printf("Writing reconsturcted map to FITS files\n");
  char outfile[100];
  sprintf(outfile, "%s%s%s%s", fileroot, "_recon_", params, ".fits");
  printf("- Outfile = %s\n", outfile);
  remove(outfile);                        // In case the file exists
  s2let_hpx_write_map(outfile, f, nside); // Now write the map to fits file

  printf("--------------------------------------------------\n");

  return 0;
}
