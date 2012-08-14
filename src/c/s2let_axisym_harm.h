// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_HARM
#define S2LET_AXISYM_HARM

void s2let_axisym_allocate_f_wav_lm(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L, int J_min);
void s2let_axisym_allocate_f_wav_multires_lm(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L, int J_min);

void s2let_axisym_allocate_wav_lm(double **wav_lm, double **scal_lm, int B, int L);
void s2let_axisym_wav_lm(double *wav_lm, double *scal_lm, int B, int L, int J_min);

void s2let_axisym_wav_analysis_lm(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);
void s2let_axisym_wav_synthesis_lm(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);

void s2let_axisym_wav_analysis_multires_lm(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);
void s2let_axisym_wav_synthesis_multires_lm(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);


#endif