// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_transform_AXISYM_LM
#define S2LET_transform_AXISYM_LM

#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

void s2let_transform_axisym_lm_allocate_f_wav(
    S2LET_COMPLEX(double) * *f_wav_lm,
    S2LET_COMPLEX(double) * *f_scal_lm,
    const s2let_parameters_t* parameters);

void s2let_transform_axisym_lm_allocate_f_wav_multires(
    S2LET_COMPLEX(double) * *f_wav_lm,
    S2LET_COMPLEX(double) * *f_scal_lm,
    const s2let_parameters_t* parameters);

void s2let_transform_axisym_lm_allocate_wav(double** wav_lm, double** scal_lm, const s2let_parameters_t* parameters);
void s2let_transform_axisym_lm_wav(double* wav_lm, double* scal_lm, const s2let_parameters_t* parameters);

void s2let_transform_axisym_lm_wav_analysis(
    S2LET_COMPLEX(double) * f_wav_lm,
    S2LET_COMPLEX(double) * f_scal_lm,
    const S2LET_COMPLEX(double) * flm,
    const double* wav_lm,
    const double* scal_lm,
    const s2let_parameters_t* parameters);
void s2let_transform_axisym_lm_wav_synthesis(
    S2LET_COMPLEX(double) * flm,
    const S2LET_COMPLEX(double) * f_wav_lm,
    const S2LET_COMPLEX(double) * f_scal_lm,
    const double* wav_lm,
    const double* scal_lm,
    const s2let_parameters_t* parameters);

void s2let_transform_axisym_lm_wav_analysis_multires(
    S2LET_COMPLEX(double) * f_wav_lm,
    S2LET_COMPLEX(double) * f_scal_lm,
    const S2LET_COMPLEX(double) * flm,
    const double* wav_lm,
    const double* scal_lm,
    const s2let_parameters_t* parameters);
void s2let_transform_axisym_lm_wav_synthesis_multires(
    S2LET_COMPLEX(double) * flm,
    const S2LET_COMPLEX(double) * f_wav_lm,
    const S2LET_COMPLEX(double) * f_scal_lm,
    const double* wav_lm,
    const double* scal_lm,
    const s2let_parameters_t* parameters);

#ifdef __cplusplus
}
#endif
#endif
