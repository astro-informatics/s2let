
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
np.import_array()

#----------------------------------------------------------------------------------------------------#

cdef extern from "ssht.h":

	double ssht_sampling_mw_t2theta(int t, int L);
	double ssht_sampling_mw_p2phi(int p, int L);
	int ssht_sampling_mw_n(int L);
	int ssht_sampling_mw_ntheta(int L);
	int ssht_sampling_mw_nphi(int L);

#----------------------------------------------------------------------------------------------------#

cdef extern from "s2let.h":

	int s2let_n_scal(const s2let_parameters_t *parameters);
	int s2let_n_wav(const s2let_parameters_t *parameters);

	int s2let_n_wav_j(int j, const s2let_parameters_t *parameters);

	void s2let_mw_alm2map(double complex * f, const double complex * flm, int L, int spin);

	void s2let_mw_map2alm(double complex * flm, const double complex * f, int L, int spin);

	int s2let_bandlimit(int j, const s2let_parameters_t *parameters);

	void s2let_synthesis_wav2lm(
		double complex *flm,
		const double complex *f_wav,
		const double complex *f_scal,
		const s2let_parameters_t *parameters
	);

	void s2let_analysis_lm2wav(
		double complex *f_wav,
		double complex *f_scal,
		const double complex *flm,
		const s2let_parameters_t *parameters
	);

	void s2let_transform_axisym_lm_allocate_wav(
		double **wav_lm, double **scal_lm, const s2let_parameters_t *parameters);

	void s2let_transform_axisym_lm_wav(
		double *wav_lm, double *scal_lm, const s2let_parameters_t *parameters);

	void s2let_transform_axisym_lm_wav_analysis(
		double complex *f_wav_lm,
		double complex *f_scal_lm,
		const double complex *flm,
		const double *wav_lm,
		const double *scal_lm,
		const s2let_parameters_t *parameters
	);

	void s2let_transform_axisym_lm_wav_synthesis(
		double complex *flm,
		const double complex *f_wav_lm,
		const double complex *f_scal_lm,
		const double *wav_lm,
		const double *scal_lm,
		const s2let_parameters_t *parameters
	);

	int s2let_j_max(s2let_parameters_t *parameters);

	ctypedef enum ssht_dl_method_t:
		SSHT_DL_RISBO, SSHT_DL_TRAPANI
	ctypedef enum s2let_sampling_t:
		S2LET_SAMPLING_MW
	ctypedef enum s2let_wav_norm_t:
		S2LET_WAV_NORM_DEFAULT, S2LET_WAV_NORM_SPIN_LOWERED

	ctypedef struct s2let_parameters_t:
		int J_min
		int B
		int L
		int N
		int upsample
		int spin
		ssht_dl_method_t dl_method
		s2let_wav_norm_t normalization;
		s2let_sampling_t sampling_scheme;
		int original_spin
		int reality
		int verbosity


#----------------------------------------------------------------------------------------------------#

cdef extern from "stdlib.h":
	void free(void* ptr)

#----------------------------------------------------------------------------------------------------#

def pys2let_j_max(B, L, J_min):
	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	return s2let_j_max(&parameters);

#----------------------------------------------------------------------------------------------------#

def healpy_lm(el, em, L):
	return em*(2*L-1-em)/2+el

#----------------------------------------------------------------------------------------------------#

def analysis_axisym_lm_wav(
	np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, B, L, J_min):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	J = s2let_j_max(&parameters);

	cdef double *wav_lm, *scal_lm;
	s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
	s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);

	f_scal_lm = np.empty([L * L,], dtype=complex)
	f_wav_lm = np.empty([L * L * (J+1-J_min),], dtype=complex)
	f_lm = lm_hp2lm(flm_hp, L)

	s2let_transform_axisym_lm_wav_analysis(
		<double complex*> np.PyArray_DATA(f_wav_lm),
		<double complex*> np.PyArray_DATA(f_scal_lm),
		<double complex*> np.PyArray_DATA(f_lm),
		wav_lm, scal_lm, &parameters
	);

	free(scal_lm);
	free(wav_lm);

	f_scal_lm_hp = np.empty([L*(L+1)/2,], dtype=complex)
	f_wav_lm_hp = np.empty([L*(L+1)/2, J-J_min+1], dtype=complex)

	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_scal_lm_hp[healpy_lm(el, em, L)] = f_scal_lm[ el * el + el + em ]
			for j from 0 <= j <= J-J_min:
				f_wav_lm_hp[healpy_lm(el, em, L), j] = f_wav_lm[ L*L*j + el * el + el + em ]

	return f_wav_lm_hp, f_scal_lm_hp

#----------------------------------------------------------------------------------------------------#

def synthesis_axisym_lm_wav(
	np.ndarray[double complex, ndim=2, mode="c"] f_wav_lm_hp not None,
	np.ndarray[double complex, ndim=1, mode="c"] f_scal_lm_hp not None, B, L, J_min):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	J = s2let_j_max(&parameters);

	cdef double *wav_lm, *scal_lm;
	s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
	s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);

	f_lm = np.empty([L * L,], dtype=complex)
	f_scal_lm = np.empty([L * L,], dtype=complex)
	f_wav_lm = np.empty([L * L * (J+1-J_min),], dtype=complex)

	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_scal_lm[ el * el + el - em ] = pow(-1.0, -em) * f_scal_lm_hp[healpy_lm(el, em, L)].conjugate()
			f_scal_lm[ el * el + el + em ] = f_scal_lm_hp[healpy_lm(el, em, L)]
			for j from 0 <= j <= J-J_min:
				f_wav_lm[ L*L*j + el * el + el - em ] = pow(-1.0, -em) * f_wav_lm_hp[healpy_lm(el, em, L), j].conjugate()
				f_wav_lm[ L*L*j + el * el + el + em ] = f_wav_lm_hp[healpy_lm(el, em, L), j]

	s2let_transform_axisym_lm_wav_synthesis(
		<double complex*> np.PyArray_DATA(f_lm),
		<double complex*> np.PyArray_DATA(f_wav_lm),
		<double complex*> np.PyArray_DATA(f_scal_lm),
		wav_lm, scal_lm, &parameters
	);

	free(scal_lm);
	free(wav_lm);
	f_lm_hp = lm2lm_hp(f_lm, L)

	return f_lm_hp

#----------------------------------------------------------------------------------------------------#

def analysis_lm2wav(
		np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None,
		B, L, J_min, N, spin, upsample):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	parameters.N = N;
	parameters.spin = spin;
	parameters.upsample = upsample;
	parameters.sampling_scheme = S2LET_SAMPLING_MW
	parameters.normalization = S2LET_WAV_NORM_DEFAULT
	parameters.dl_method = SSHT_DL_RISBO
	parameters.original_spin = 0
	parameters.reality = 0
	parameters.verbosity = 0

	f_scal = np.empty([s2let_n_scal(&parameters),], dtype=complex)
	f_wav = np.empty([s2let_n_wav(&parameters),], dtype=complex)
	f_lm = lm_hp2lm(flm_hp, L)

	s2let_analysis_lm2wav(
		<double complex*> np.PyArray_DATA(f_wav),
		<double complex*> np.PyArray_DATA(f_scal),
		<const double complex*> np.PyArray_DATA(f_lm),
		&parameters);

	return f_wav, f_scal

#----------------------------------------------------------------------------------------------------#

def synthesis_wav2lm(
		np.ndarray[double complex, ndim=1, mode="c"] f_wav not None,
		np.ndarray[double complex, ndim=1, mode="c"] f_scal not None,
		B, L, J_min, N, spin, upsample):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	parameters.N = N;
	parameters.spin = spin;
	parameters.upsample = upsample;
	parameters.sampling_scheme = S2LET_SAMPLING_MW
	parameters.normalization = S2LET_WAV_NORM_DEFAULT
	parameters.dl_method = SSHT_DL_RISBO
	parameters.original_spin = 0
	parameters.reality = 0
	parameters.verbosity = 0

	f_lm = np.empty([L * L,], dtype=complex)
	s2let_synthesis_wav2lm(
		<double complex*> np.PyArray_DATA(f_lm),
		<const double complex*> np.PyArray_DATA(f_wav),
		<const double complex*> np.PyArray_DATA(f_scal),
		&parameters);
	f_lm_hp = lm2lm_hp(f_lm, L)

	return f_lm_hp

#----------------------------------------------------------------------------------------------------#

def mw_size(L):
	return L*(2*L-1)

#----------------------------------------------------------------------------------------------------#

def lm2lm_hp(np.ndarray[double complex, ndim=1, mode="c"] f_lm not None, L):

	f_lm_hp = np.empty([L*(L+1)/2,], dtype=complex)
	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_lm_hp[healpy_lm(el, em, L)] = f_lm[ el * el + el + em ]

	return f_lm_hp

#----------------------------------------------------------------------------------------------------#

def lm_hp2lm(np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, L):

	f_lm = np.empty([L*L,], dtype=complex)
	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_lm[ el * el + el - em ] = pow(-1.0, -em) * ( flm_hp[healpy_lm(el, em, L)] ).conjugate()
			f_lm[ el * el + el + em ] = flm_hp[healpy_lm(el, em, L)]

	return f_lm

#----------------------------------------------------------------------------------------------------#

def wav_ind(j, n, B, L, N, J_min, upsample):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	parameters.N = N;
	parameters.spin = 0;
	parameters.upsample = upsample;
	parameters.sampling_scheme = S2LET_SAMPLING_MW
	parameters.normalization = S2LET_WAV_NORM_DEFAULT
	parameters.dl_method = SSHT_DL_RISBO
	parameters.original_spin = 0
	parameters.reality = 0
	parameters.verbosity = 0

	offset = 0
	for jprime from J_min <= jprime < j:
		offset += s2let_n_wav_j(jprime, &parameters);
	if upsample:
		bandlimit = L
	else:
		bandlimit = min(s2let_bandlimit(j, &parameters), L)
	for en from 0 <= en < n:
		offset += mw_size(bandlimit)
	nelem_wav = s2let_n_wav_j(j, &parameters);

	return offset, bandlimit, mw_size(bandlimit), nelem_wav

#----------------------------------------------------------------------------------------------------#

def alm2map_mw(np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, L, spin):

	f_lm = lm_hp2lm(flm_hp, L)
	f = np.empty([mw_size(L),], dtype=complex)
	s2let_mw_alm2map(
		<double complex*> np.PyArray_DATA(f),
		<const double complex*> np.PyArray_DATA(f_lm),
		L, spin)

	return f

#----------------------------------------------------------------------------------------------------#

def map2alm_mw(np.ndarray[double complex, ndim=1, mode="c"] f not None, L, spin):

	f_lm = np.empty([L * L,], dtype=complex)
	s2let_mw_map2alm(
		<double complex*> np.PyArray_DATA(f_lm),
		<const double complex*> np.PyArray_DATA(f),
		L, spin)
	f_lm_hp = lm2lm_hp(f_lm, L)

	return f_lm_hp

#----------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------#

def mw_sampling(L):

	ntheta = ssht_sampling_mw_ntheta(L)
	nphi = ssht_sampling_mw_nphi(L)
	thetas = np.empty([ntheta,], dtype=float)
	phis = np.empty([nphi,], dtype=float)
	for t from 0<=t<ntheta:
		thetas[t] = ssht_sampling_mw_t2theta(t, L)
	for p from 0<=p<nphi:
		phis[p] = ssht_sampling_mw_p2phi(p, L)

	return thetas, phis


#----------------------------------------------------------------------------------------------------#

