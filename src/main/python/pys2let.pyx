
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
np.import_array()

#----------------------------------------------------------------------------------------------------#

cdef extern from "so3.h":

	int so3_sampling_f_size(so3_parameters_t *params);

	ctypedef struct so3_parameters_t:
		int reality
		int L0
		int L
		int N
		so3_sampling_t sampling_scheme;
		ssht_dl_method_t dl_method;

	ctypedef enum so3_sampling_t:
		SO3_SAMPLING_MW, SO3_SAMPLING_MW_SS, SO3_SAMPLING_SIZE

#----------------------------------------------------------------------------------------------------#

cdef extern from "ssht.h":

	double ssht_sampling_mw_t2theta(int t, int L);
	double ssht_sampling_mw_p2phi(int p, int L);
	int ssht_sampling_mw_n(int L);
	int ssht_sampling_mw_ntheta(int L);
	int ssht_sampling_mw_nphi(int L);
	double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size);
	ctypedef enum ssht_dl_size_t:
		SSHT_DL_QUARTER_EXTENDED, SSHT_DL_HALF, SSHT_DL_FULL
	void ssht_dl_beta_risbo_full_table(double *dl, double beta, int L,
					   ssht_dl_size_t dl_size,
					   int el, double *sqrt_tbl);

	void ssht_dl_beta_risbo_half_table(double *dl, double beta, int L,
					   ssht_dl_size_t dl_size,
					   int el, double *sqrt_tbl, double *signs);

#----------------------------------------------------------------------------------------------------#

cdef extern from "s2let.h":

	void fill_so3_parameters(so3_parameters_t *parameters1, s2let_parameters_t *parameters2)
	int s2let_n_phi(const s2let_parameters_t *parameters);
	int s2let_n_theta(const s2let_parameters_t *parameters);
	void s2let_tiling_direction_allocate(double complex **s_elm, const s2let_parameters_t *parameters);
	void s2let_tiling_direction(double complex *s_elm, const s2let_parameters_t *parameters);

	void s2let_tiling_wavelet_allocate(double complex **psi, double **phi, const s2let_parameters_t *parameters);
	void s2let_tiling_wavelet(double complex *psi, double *phi, const s2let_parameters_t *parameters);
	void s2let_tiling_axisym(double *kappa, double *kappa0, const s2let_parameters_t *parameters);
	void s2let_tiling_axisym_allocate(double **kappa, double **kappa0, const s2let_parameters_t *parameters);

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

	void s2let_analysis_lm2wav_manual(
		double complex *f_wav,
		double complex *f_scal,
		const double complex *flm,
		const double *scal_l,
		const double complex *wav_lm,
		const int scal_bandlimit,
		const int *wav_bandlimits,
		int J,
		int L,
		int spin,
		int N
	);

	void s2let_synthesis_wav2lm_manual(
		double complex *flm,
		const double complex *f_wav,
		const double complex *f_scal,
		const double *scal_l,
		const double complex *wav_lm,
		const int scal_bandlimit,
		const int *wav_bandlimits,
		int J,
		int L,
		int spin,
		int N
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

	ctypedef struct s2let_parameters_t:
		int J_min
		double B
		int L
		int N
		int upsample
		int spin
		ssht_dl_method_t dl_method
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

def mw_lm(el, em):
	return el*el + el * em

#----------------------------------------------------------------------------------------------------#

def healpy_lm(el, em, L):
	return em*(2*L-1-em)/2+el

#----------------------------------------------------------------------------------------------------#

def analysis_axisym_lm_wav(
	np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, B, L, J_min, spin_lowered = False):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	J = s2let_j_max(&parameters);

	cdef double *wav_lm, *scal_lm;
	s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
	s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);

	f_scal_lm = np.zeros([L * L,], dtype=complex)
	f_wav_lm = np.zeros([L * L * (J+1-J_min),], dtype=complex)
	f_lm = lm_hp2lm(flm_hp, L)

	s2let_transform_axisym_lm_wav_analysis(
		<double complex*> np.PyArray_DATA(f_wav_lm),
		<double complex*> np.PyArray_DATA(f_scal_lm),
		<double complex*> np.PyArray_DATA(f_lm),
		wav_lm, scal_lm, &parameters
	);

	free(scal_lm);
	free(wav_lm);

	f_scal_lm_hp = np.zeros([L*(L+1)/2,], dtype=complex)
	f_wav_lm_hp = np.zeros([L*(L+1)/2, J-J_min+1], dtype=complex)

	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_scal_lm_hp[healpy_lm(el, em, L)] = f_scal_lm[ el * el + el + em ]
			for j from 0 <= j <= J-J_min:
				f_wav_lm_hp[healpy_lm(el, em, L), j] = f_wav_lm[ L*L*j + el * el + el + em ]

	return f_wav_lm_hp, f_scal_lm_hp

#----------------------------------------------------------------------------------------------------#

def synthesis_axisym_lm_wav(
	np.ndarray[double complex, ndim=2, mode="c"] f_wav_lm_hp not None,
	np.ndarray[double complex, ndim=1, mode="c"] f_scal_lm_hp not None, B, L, J_min, spin_lowered = False):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	J = s2let_j_max(&parameters);

	cdef double *wav_lm, *scal_lm;
	s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
	s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);

	f_lm = np.zeros([L * L,], dtype=complex)
	f_scal_lm = np.zeros([L * L,], dtype=complex)
	f_wav_lm = np.zeros([L * L * (J+1-J_min),], dtype=complex)

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

def verify_tiling(L,
		np.ndarray[double, ndim=1, mode="c"] scal_l,
		np.ndarray[double, ndim=1, mode="c"] wav_l,
		scal_bandlimit,
		np.ndarray[int, ndim=1, mode="c"] wav_bandlimits):

	J = wav_bandlimits.size - 1

	totalsum = np.zeros((L,))

	totalsum += (scal_l)**2.0
	cumsum = np.cumsum(scal_l)
	cumsum /= cumsum[-1]
	if cumsum[scal_bandlimit-1] < 1:
		print 'Wrong band-limit given for scaling function (', scal_bandlimit, ')'
		print 'Measured band-limit is ', np.where(cumsum == 1)[0][0]+1
		return False
	print 'Measured / given band-limits for scaling fct are', np.where(cumsum == 1)[0][0]+1, '/', scal_bandlimit

	ell = np.arange(L)
	for j from 0 <= j <= J:
		totalsum += (wav_l[j*L:(j+1)*L])**2.0
		cumsum = np.cumsum(wav_l[j*L:(j+1)*L])
		cumsum /= cumsum[-1]
		if cumsum[wav_bandlimits[j]-1] < 1:
			print 'Wrong band-limit given for wavelet', j+1, 'on', J+1, '(', wav_bandlimits[j], ')'
			print 'Measured band-limit is ', np.where(cumsum == 1)[0][0]+1
			print 'Wavelet:',  wav_l[j*L:(j+1)*L]
			return False
		print 'Measured / given band-limits for wavelet', j, 'are', np.where(cumsum == 1)[0][0]+1, '/', wav_bandlimits[j]

	if np.allclose(totalsum, 1.0):
		return True
	else:
		print 'Admissibility condition not satisfied:'
		print totalsum
		return False

#----------------------------------------------------------------------------------------------------#

def analysis_lm2wav_manualtiling(
		np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None,
		L, N, spin,
		np.ndarray[double, ndim=1, mode="c"] kappa_scal_l,
		np.ndarray[double, ndim=1, mode="c"] kappa_wav_l,
		scal_bandlimit,
		np.ndarray[int, ndim=1, mode="c"] wav_bandlimits):

	J = wav_bandlimits.size - 1

	cdef s2let_parameters_t parameters = {};
	parameters.L = L;
	parameters.N = N;
	parameters.spin = spin;
	parameters.upsample = 0;
	parameters.J_min = 0;

	cdef so3_parameters_t so3_parameters = {};
	fill_so3_parameters(&so3_parameters, &parameters);
	so3_parameters.sampling_scheme = SO3_SAMPLING_MW;

	scal_l = np.zeros([L,], dtype=float)
	for el from 0 <= el < L:
		scal_l[el] = np.sqrt((2*el+1)/(4.0*np.pi)) * kappa_scal_l[el];

	wav_lm = np.zeros([L*L*(J+1),], dtype=complex)
	cdef double complex *s_elm;
	s2let_tiling_direction_allocate(&s_elm, &parameters);
	s2let_tiling_direction(s_elm, &parameters);
	for j from 0 <= j <= J:
		for el from 0 <= el < L:
			for em from -el <= em <= el:
				ind = el*el + el + em
				wav_lm[j*L*L + ind] = np.sqrt((2*el+1)/(8.0*np.pi*np.pi)) * kappa_wav_l[j*L + el] * s_elm[ind];
	free(s_elm);

	cdef s2let_parameters_t bl_parameters = {};
	bl_parameters.L = int(scal_bandlimit);
	total = s2let_n_phi(&bl_parameters) * s2let_n_theta(&bl_parameters)
	f_scal = np.zeros([total,], dtype=complex)
	total = 0
	for j from 0 <= j <= J:
		so3_parameters.L0 = 0;
		bandlimit = np.min([wav_bandlimits[j], L]);
		so3_parameters.L = bandlimit;
		Nj = np.min([N, bandlimit]);
		Nj += (Nj+N) % 2; #// ensure N and Nj are both even or both odd
		so3_parameters.N = Nj;
		total += so3_sampling_f_size(&so3_parameters)
	f_wav = np.zeros([total,], dtype=complex)

	f_lm = lm_hp2lm(flm_hp, L)

	print 'scal_bandlimit = ', int(scal_bandlimit)
	print 'wav_bandlimits = ', wav_bandlimits
	print 'J, L, spin, N = ', J, L, spin, N

	print 'Done pre-computing - running s2let_analysis_lm2wav_manual'
	s2let_analysis_lm2wav_manual(
		<double complex*> np.PyArray_DATA(f_wav),
		<double complex*> np.PyArray_DATA(f_scal),
		<const double complex*> np.PyArray_DATA(f_lm),
		<const double*> np.PyArray_DATA(scal_l),
		<const double complex*> np.PyArray_DATA(wav_lm),
		scal_bandlimit,
		<const int*> np.PyArray_DATA(wav_bandlimits),
		J, L, spin, N);

	return f_wav, f_scal

#----------------------------------------------------------------------------------------------------#

def synthesis_wav2lm_manualtiling(
		np.ndarray[double complex, ndim=1, mode="c"] f_wav not None,
		np.ndarray[double complex, ndim=1, mode="c"] f_scal not None,
		L, N, spin,
		np.ndarray[double, ndim=1, mode="c"] kappa_scal_l,
		np.ndarray[double, ndim=1, mode="c"] kappa_wav_l,
		scal_bandlimit,
		np.ndarray[int, ndim=1, mode="c"] wav_bandlimits):

	J = wav_bandlimits.size - 1

	cdef s2let_parameters_t parameters = {};
	parameters.L = L;
	parameters.N = N;
	parameters.spin = spin;
	parameters.upsample = 0;
	parameters.J_min = 0;

	scal_l = np.zeros([L,], dtype=float)
	for el from 0 <= el < L:
		scal_l[el] = np.sqrt((2*el+1)/(4.0*np.pi)) * kappa_scal_l[el];

	wav_lm = np.zeros([L*L*(J+1),], dtype=complex)
	cdef double complex *s_elm;
	s2let_tiling_direction_allocate(&s_elm, &parameters);
	s2let_tiling_direction(s_elm, &parameters);
	for j from 0 <= j <= J:
		for el from 0 <= el < L:
			for em from -el <= em <= el:
				ind = el*el + el + em
				wav_lm[j*L*L + ind] = np.sqrt((2*el+1)/(8.0*np.pi*np.pi)) * kappa_wav_l[j*L + el] * s_elm[ind];
	free(s_elm);

	f_lm = np.zeros([L * L,], dtype=complex)

	print 'Done pre-computing - running s2let_synthesis_wav2lm_manual'
	s2let_synthesis_wav2lm_manual(
		<double complex*> np.PyArray_DATA(f_lm),
		<const double complex*> np.PyArray_DATA(f_wav),
		<const double complex*> np.PyArray_DATA(f_scal),
		<const double*> np.PyArray_DATA(scal_l),
		<const double complex*> np.PyArray_DATA(wav_lm),
		scal_bandlimit,
		<const int*> np.PyArray_DATA(wav_bandlimits),
		J, L, spin, N);

	f_lm_hp = lm2lm_hp(f_lm, L)

	return f_lm_hp

#----------------------------------------------------------------------------------------------------#

def analysis_lm2wav(
		np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None,
		B, L, J_min, N, spin, upsample, spin_lowered=False, original_spin=0):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	parameters.N = N;
	parameters.spin = spin;
	parameters.upsample = upsample;
	parameters.sampling_scheme = S2LET_SAMPLING_MW
	parameters.original_spin = original_spin
	parameters.dl_method = SSHT_DL_RISBO
	parameters.reality = 0
	parameters.verbosity = 0

	f_scal = np.zeros([s2let_n_scal(&parameters),], dtype=complex)
	f_wav = np.zeros([s2let_n_wav(&parameters),], dtype=complex)
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
		B, L, J_min, N, spin, upsample, spin_lowered=False, original_spin=0):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	parameters.N = N;
	parameters.spin = spin;
	parameters.upsample = upsample;
	parameters.sampling_scheme = S2LET_SAMPLING_MW
	parameters.original_spin = original_spin
	parameters.dl_method = SSHT_DL_RISBO
	parameters.reality = 0
	parameters.verbosity = 0

	f_lm = np.zeros([L * L,], dtype=complex)
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

	f_lm_hp = np.zeros([L*(L+1)/2,], dtype=complex)
	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_lm_hp[healpy_lm(el, em, L)] = f_lm[ el * el + el + em ]

	return f_lm_hp

#----------------------------------------------------------------------------------------------------#

def lm_hp2lm(np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, L):

	f_lm = np.zeros([L*L,], dtype=complex)
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

def mw_sampling(L):

	ntheta = ssht_sampling_mw_ntheta(L)
	nphi = ssht_sampling_mw_nphi(L)
	thetas = np.zeros([ntheta,], dtype=float)
	phis = np.zeros([nphi,], dtype=float)
	for t from 0<=t<ntheta:
		thetas[t] = ssht_sampling_mw_t2theta(t, L)
	for p from 0<=p<nphi:
		phis[p] = ssht_sampling_mw_p2phi(p, L)

	return thetas, phis


#----------------------------------------------------------------------------------------------------#

def alm2map_mw(np.ndarray[double complex, ndim=1, mode="c"] f_lm not None, L, spin):

	#f_lm = lm_hp2lm(f_lm_hp, L)
	f = np.zeros([mw_size(L),], dtype=complex)
	s2let_mw_alm2map(
		<double complex*> np.PyArray_DATA(f),
		<const double complex*> np.PyArray_DATA(f_lm),
		L, spin)

	return f

#----------------------------------------------------------------------------------------------------#

def map2alm_mw(np.ndarray[double complex, ndim=1, mode="c"] f not None, L, spin):

	f_lm = np.zeros([L * L,], dtype=complex)
	s2let_mw_map2alm(
		<double complex*> np.PyArray_DATA(f_lm),
		<const double complex*> np.PyArray_DATA(f),
		L, spin)
	#f_lm_hp = lm2lm_hp(f_lm, L)

	return f_lm

#----------------------------------------------------------------------------------------------------#

def wavelet_tiling(B, L, N, J_min, spin, original_spin=0):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.N = N;
	parameters.spin = spin;
	parameters.original_spin = 0
	parameters.J_min = J_min;
	J = s2let_j_max(&parameters);

	cdef double complex *psi;
	cdef double *phi;
	s2let_tiling_wavelet_allocate(&psi, &phi, &parameters)
	s2let_tiling_wavelet(psi, phi, &parameters);

	scal_l = np.zeros([L,], dtype=complex)
	wav_l = np.zeros([L*L, J-J_min+1], dtype=complex)

	for el from 0 <= el < L:
		scal_l[el] = phi[el];

	for el from 0 <= el < L*L:
		for j from J_min <= j <= J:
			wav_l[el, j-J_min] = psi[ j*L*L + el];

	free(psi);
	free(phi);

	return scal_l, wav_l

#----------------------------------------------------------------------------------------------------#

def axisym_wav_l(B, L, J_min):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	J = s2let_j_max(&parameters);

	cdef double *kappa, *kappa0;
	s2let_tiling_axisym_allocate(&kappa, &kappa0, &parameters);
	s2let_tiling_axisym(kappa, kappa0, &parameters);

	scal_l = np.zeros([L,], dtype=float)
	wav_l = np.zeros([L, J-J_min+1], dtype=float)

	for el from 0 <= el < L:
		scal_l[el] = kappa0[el];
		for j from 0 <= j <= J-J_min:
			wav_l[el, j] = kappa[el+(j+J_min)*L];

	free(kappa);
	free(kappa0);

	return scal_l, wav_l

#----------------------------------------------------------------------------------------------------#

def ssht_dl_beta_risbo(beta, L):

	dl_beta = np.zeros([L,(2*L-1),(2*L-1)], dtype=float)
	cdef ssht_dl_size_t ssht_dl_size = SSHT_DL_FULL

	#sqrt_tbl = np.zeros([2*(L-1)+2,], dtype=float)
	#for el from 0 <= el < 2*(L-1)+1:
	#	sqrt_tbl[el] = np.sqrt(el);
	sqrt_tbl = np.sqrt(np.arange(2.0*(L-1)+2.0))

	#signs = np.ones([L+1,], dtype=float)
	signs = np.power(-1.0, np.arange(L+1));

	#cdef double *dl;
	#dl = ssht_dl_calloc(L, SSHT_DL_FULL);
	dl = np.zeros([(2*L-1)*(2*L-1)], dtype=float)

	#dl_offset = L - 1; #//ssht_dl_get_offset(L, SSHT_DL_FULL);
	#dl_stride = 2*L - 1; #//ssht_dl_get_stride(L, SSHT_DL_FULL);

	for el from 0 <= el < L:
		ssht_dl_beta_risbo_half_table(
			<double*> np.PyArray_DATA(dl),
			beta, L, SSHT_DL_FULL, el,
			<double*> np.PyArray_DATA(sqrt_tbl),
			<double*> np.PyArray_DATA(signs));
		dl_beta[el,:,:] = dl.reshape((2*L-1, 2*L-1))

	#free(dl);

	return dl_beta

#----------------------------------------------------------------------------------------------------#

# Function to construct a hybrid wavelet tiling that should be valid for an invertible wavelet transform
def construct_hybrid_tiling(L, J_min0, Bs, L_transitions):
	nb = Bs.size
	J_mins = np.repeat(0, nb)
	J_mins[0] = J_min0
	Js = [pys2let_j_max(B, L, J_min) for B, J_min in zip(Bs, J_mins)]
	J = 0
	L_bounds = np.zeros((nb+1,) , dtype=np.int32)
	j_transitions_left = np.zeros((nb,) , dtype=np.int32)
	j_transitions_right = np.zeros((nb,) , dtype=np.int32)
	for k in range(nb):
		if k == 0:
			j_transitions_left[k] = J_mins[k]
		else:
			j_transitions_left[k] = int(np.log(L_bounds[k]) / np.log(Bs[k])) + 1
		if k == nb - 1:
			j_transitions_right[k] = Js[k] + 1
		else:
			j_transitions_right[k] = int(np.log(L_transitions[k]) / np.log(Bs[k])) + 1
		if k < nb - 1:
			L_bounds[k+1] = np.rint(Bs[k]**(j_transitions_right[k] - 1))
		J += int( j_transitions_right[k] - j_transitions_left[k] - 1 )
	L_bounds[-1] = L
	hybrid_wav_l = np.zeros((L, J+1))
	hybrid_wav_bandlimits = np.zeros((J+1,), dtype=np.int32)
	off = 0
	for k in range(nb):
		#print 'k =',k+1,'on',nb,'with j_left =', j_transitions_left[k],' and j_right =', j_transitions_right[k]
		scal_l, wav_l = axisym_wav_l(Bs[k], L, J_mins[k])
		jrange = np.arange(J_mins[k],Js[k]+1)
		wav_bandlimits = np.array([np.rint(np.min([x,L])) for x in Bs[k]**(jrange+1)]).astype(np.int32)
		if k == 0:
			hybrid_scal_bandlimit = Bs[k]**(J_mins[k])
			hybrid_scal_l = np.zeros((L, ))
			hybrid_scal_l[:] = scal_l[:]

		for j in range(j_transitions_left[k]-J_mins[k],j_transitions_right[k]-J_mins[k]):
			hybrid_wav_bandlimits[off] = wav_bandlimits[j]
			hybrid_wav_l[L_bounds[k]:L_bounds[k+1],off] = wav_l[L_bounds[k]:L_bounds[k+1],j]
			if j == j_transitions_left[k]- J_mins[k] and wav_l[L_bounds[k],j] < 1.0:
				hybrid_wav_l[L_bounds[k]:L_bounds[k+1],off] = np.sqrt(hybrid_wav_l[L_bounds[k]:L_bounds[k+1],off]**2 + wav_l[L_bounds[k]:L_bounds[k+1],j-1]**2)
			if k < nb -1 and j == j_transitions_right[k] - 1 - J_mins[k] and wav_l[L_bounds[k],j] < 1.0:
				hybrid_wav_l[L_bounds[k]:L_bounds[k+1],off] = np.sqrt(hybrid_wav_l[L_bounds[k]:L_bounds[k+1],off]**2 + wav_l[L_bounds[k]:L_bounds[k+1],j+1]**2)
			if j < j_transitions_right[k]-1- J_mins[k]:
				off += 1

	return hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits, J, L_bounds

#----------------------------------------------------------------------------------------------------#
