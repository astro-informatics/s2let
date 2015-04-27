#ifndef S2LET_HELPER
#define S2LET_HELPER

int s2let_n_phi(const s2let_parameters_t *parameters);
int s2let_n_theta(const s2let_parameters_t *parameters);
int s2let_n_px(const s2let_parameters_t *parameters);

int s2let_n_lm(const s2let_parameters_t *parameters);

int s2let_n_lmn(const s2let_parameters_t *parameters);
int s2let_n_lm_scal(const s2let_parameters_t *parameters);
int s2let_n_lmn_wav(const s2let_parameters_t *parameters);

int s2let_n_gamma(const s2let_parameters_t *parameters);
int s2let_n_scal(const s2let_parameters_t *parameters);
int s2let_n_wav(const s2let_parameters_t *parameters);
int s2let_n_wav_j(int j, const s2let_parameters_t *parameters);

int s2let_n_lmn_cur(const s2let_parameters_t *parameters);
int s2let_n_cur(const s2let_parameters_t *parameters);
int s2let_n_cur_j(int j, const s2let_parameters_t *parameters);

#endif
