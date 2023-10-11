#include <ssht/ssht.h>

#include "s2let/s2let.h"

int s2let_n_phi(const s2let_parameters_t *parameters) {
  if (parameters->sampling_scheme == S2LET_SAMPLING_MW_SS)
    return ssht_sampling_mw_ss_nphi(parameters->L);
  else
    return ssht_sampling_mw_nphi(parameters->L);
}

int s2let_n_theta(const s2let_parameters_t *parameters) {
  if (parameters->sampling_scheme == S2LET_SAMPLING_MW_SS)
    return ssht_sampling_mw_ss_ntheta(parameters->L);
  else
    return ssht_sampling_mw_ntheta(parameters->L);
}

int s2let_n_px(const s2let_parameters_t *parameters) {
  return s2let_n_phi(parameters) * s2let_n_theta(parameters);
}

int s2let_n_lm(const s2let_parameters_t *parameters) {
  int L = parameters->L;
  return L * L;
}

int s2let_n_lm_scal(const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;
  int bandlimit =
      (parameters->upsample) ? L : MIN(s2let_bandlimit(J_min - 1, parameters), L);

  return bandlimit * bandlimit;
}

int s2let_n_lmn_wav(const s2let_parameters_t *parameters) {
  so3_parameters_t so3_parameters = {};
  fill_so3_parameters(&so3_parameters, parameters);

  int L = parameters->L;
  int N = parameters->N;
  int J_min = parameters->J_min;
  int J = s2let_j_max(parameters);
  int bandlimit = L;
  int j, total = 0;
  for (j = J_min; j <= J; ++j) {
    if (!parameters->upsample) {
      bandlimit = MIN(s2let_bandlimit(j, parameters), L);
      so3_parameters.L = bandlimit;
      so3_parameters.N = MIN(N, bandlimit);
    }
    total += so3_sampling_flmn_size(&so3_parameters);
  }
  return total;
}

int s2let_n_gamma(const s2let_parameters_t *parameters) {
  so3_parameters_t so3_parameters = {};
  fill_so3_parameters(&so3_parameters, parameters);

  return so3_sampling_ngamma(&so3_parameters);
}

int s2let_n_scal(const s2let_parameters_t *parameters) {
  int J_min = parameters->J_min;
  int L = parameters->L;
  int bandlimit = (parameters->upsample)
                      ? parameters->L
                      : MIN(s2let_bandlimit(J_min - 1, parameters), L);

  s2let_parameters_t bl_parameters = *parameters;
  bl_parameters.L = bandlimit;

  return s2let_n_phi(&bl_parameters) * s2let_n_theta(&bl_parameters);
}

int s2let_n_wav(const s2let_parameters_t *parameters) {
  so3_parameters_t so3_parameters = {};
  fill_so3_parameters(&so3_parameters, parameters);

  int L = parameters->L;
  int J_min = parameters->J_min;
  int J = s2let_j_max(parameters);
  int bandlimit = L;
  int j, total = 0;
  for (j = J_min; j <= J; ++j) {
    if (!parameters->upsample) {
      bandlimit = MIN(s2let_bandlimit(j, parameters), L);
      so3_parameters.L = bandlimit;
    }
    total += so3_sampling_f_size(&so3_parameters);
  }
  return total;
}

int s2let_n_wav_j(int j, const s2let_parameters_t *parameters) {
  so3_parameters_t so3_parameters = {};
  fill_so3_parameters(&so3_parameters, parameters);

  if (!parameters->upsample) {
    int L = parameters->L;
    so3_parameters.L = MIN(s2let_bandlimit(j, parameters), L);
  }

  return so3_sampling_f_size(&so3_parameters);
}
