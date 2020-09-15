// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_SYNTHESIS
#define S2LET_SYNTHESIS

#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Harmonic-space wavelet transform **/

void s2let_synthesis_lmn2lm(
    S2LET_COMPLEX(double) * flm,
    const S2LET_COMPLEX(double) * f_wav_lmn,
    const S2LET_COMPLEX(double) * f_scal_lm,
    const S2LET_COMPLEX(double) * wav_lm,
    const double* scal_l,
    const s2let_parameters_t* parameters);

/** Harmonic-space wavelet transform for real signals **/

void s2let_synthesis_lmn2lm_real(
    S2LET_COMPLEX(double) * flm,
    const S2LET_COMPLEX(double) * f_wav_lmn,
    const S2LET_COMPLEX(double) * f_scal_lm,
    const S2LET_COMPLEX(double) * wav_lm,
    const double* scal_l,
    const s2let_parameters_t* parameters);

/** Harmonic-space wavelet transform, manual tiling **/

void s2let_synthesis_wav2lm_manual(
    S2LET_COMPLEX(double) * flm,
    const S2LET_COMPLEX(double) * f_wav,
    const S2LET_COMPLEX(double) * f_scal,
    const double* scal_l,
    const S2LET_COMPLEX(double) * wav_lm,
    const int scal_bandlimit,
    const int* wav_bandlimits,
    int J,
    int L,
    int spin,
    int N);

void s2let_synthesis_wav2lm(
    S2LET_COMPLEX(double) * flm,
    const S2LET_COMPLEX(double) * f_wav,
    const S2LET_COMPLEX(double) * f_scal,
    const s2let_parameters_t* parameters);

/** Harmonic-space wavelet transform for real signals **/

void s2let_synthesis_wav2lm_real(
    S2LET_COMPLEX(double) * flm,
    const double* f_wav,
    const double* f_scal,
    const s2let_parameters_t* parameters);

/** Pixel-space wavelet transform **/

void s2let_synthesis_wav2px(
    S2LET_COMPLEX(double) * f,
    const S2LET_COMPLEX(double) * f_wav,
    const S2LET_COMPLEX(double) * f_scal,
    const s2let_parameters_t* parameters);

/** Pixel-space wavelet transform for real signals**/

void s2let_synthesis_wav2px_real(
    double* f,
    const double* f_wav,
    const double* f_scal,
    const s2let_parameters_t* parameters);

#ifdef __cplusplus
}
#endif
#endif
