// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_ANALYSIS
#define S2LET_ANALYSIS

#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Harmonic-space wavelet transform **/

void s2let_analysis_lm2lmn(
    S2LET_COMPLEX(double) * f_wav_lmn,
    S2LET_COMPLEX(double) * f_scal_lm,
    const S2LET_COMPLEX(double) * flm,
    const S2LET_COMPLEX(double) * wav_lm,
    const double* scal_l,
    const s2let_parameters_t* parameters);

/** Harmonic-space wavelet transform for real signals **/

void s2let_analysis_lm2lmn_real(
    S2LET_COMPLEX(double) * f_wav_lmn,
    S2LET_COMPLEX(double) * f_scal_lm,
    const S2LET_COMPLEX(double) * flm,
    const S2LET_COMPLEX(double) * wav_lm,
    const double* scal_l,
    const s2let_parameters_t* parameters);

/** Harmonic-space wavelet transform **/

void s2let_analysis_lm2wav_manual(
    S2LET_COMPLEX(double) * f_wav,
    S2LET_COMPLEX(double) * f_scal,
    const S2LET_COMPLEX(double) * flm,
    const double* scal_l,
    const S2LET_COMPLEX(double) * wav_lm,
    const int scal_bandlimit,
    const int* wav_bandlimits,
    int J,
    int L,
    int spin,
    int N);

void s2let_analysis_lm2wav(
    S2LET_COMPLEX(double) * f_wav,
    S2LET_COMPLEX(double) * f_scal,
    const S2LET_COMPLEX(double) * flm,
    const s2let_parameters_t* parameters);

/** Harmonic-space wavelet transform for real signals **/

void s2let_analysis_lm2wav_real(
    double* f_wav,
    double* f_scal,
    const S2LET_COMPLEX(double) * flm,
    const s2let_parameters_t* parameters);

/** Pixel-space wavelet transform **/

void s2let_analysis_px2wav(
    S2LET_COMPLEX(double) * f_wav,
    S2LET_COMPLEX(double) * f_scal,
    const S2LET_COMPLEX(double) * f,
    const s2let_parameters_t* parameters);

/** Pixel-space wavelet transform for real signals**/

void s2let_analysis_px2wav_real(
    double* f_wav,
    double* f_scal,
    const double* f,
    const s2let_parameters_t* parameters);

#ifdef __cplusplus
}
#endif
#endif
