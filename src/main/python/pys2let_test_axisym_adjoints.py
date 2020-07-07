from pys2let import *
import numpy as np

# DOT PRODUCT TESTS
# if y = Ax and g = A'f, show that f'y = g'x
L = 10
B = 2
J_min = 2
J = pys2let_j_max(B, L, J_min)
nwvlts = J - J_min + 1

def random_lms(L):
    lms = np.zeros(L * L, dtype=np.complex)
    for el in range(L):
        em = 0
        while em <= el:
            rand = np.asarray(np.random.rand(), dtype=np.complex)
            lms[el * el + el - em] = pow(-1.0, -em) * rand.conjugate()
            lms[el * el + el + em] = rand
            em += 1
    return lms

def random_mw_map(L):
    f_lm = random_lms(L)
    f = alm2map_mw(f_lm, L, 0)
    return f

def random_mw_wavelet_maps(L):
    f_wav_mw = np.column_stack([random_mw_map(L) for _ in range(nwvlts)])
    f_wav_mw = np.ascontiguousarray(f_wav_mw.flatten("F"))
    f_scal_mw = np.ascontiguousarray(random_mw_map(L))
    return f_scal_mw, f_wav_mw


#  ANALYSIS
# fwd: input = map, output = wavelet/scaling coeffs
# some random input map
x = random_mw_map(L)  # MW map

# perform fwd transform via hp
x_lm = map2alm_mw(x, L, 0)
x_hp_lm = lm2lm_hp(x_lm, L)
y_wav_hp_lm, y_scal_hp_lm = analysis_axisym_lm_wav(x_hp_lm, B, L, J_min)
y_scal_mw = alm2map_mw(lm_hp2lm(y_scal_hp_lm, L), L, 0)
y_wav_mw = np.zeros((mw_size(L) * nwvlts), dtype=np.complex)
offset = mw_size(L)
for i in range(nwvlts):
    y_wav_mw[i*offset:(i+1)*offset] = alm2map_mw(
        lm_hp2lm(np.ascontiguousarray(y_wav_hp_lm[:,i]), L), 
        L, 0)
y = np.concatenate((y_scal_mw, y_wav_mw))  # MW 

# adj: input = wavelet/scaling coeffs, output = map
# some random wavelet coeffs
f_scal_mw, f_wav_mw = random_mw_wavelet_maps(L)
f = np.concatenate((f_scal_mw, f_wav_mw))  # MW 

# perform adjoint transform
g = analysis_adjoint_axisym_wav_mw(f_wav_mw, f_scal_mw, B, L, J_min)  # MW map

LHS = f.T.dot(y)
RHS = g.T.dot(x)
print(f"Analysis: {LHS}, {RHS}")
print(f"Dot product test error: {LHS - RHS}")


# SYNTHESIS
# inv: input = wavelet/scaling coeffs, output = map 
# some random wavelet coeffs
x_scal_mw, x_wav_mw = random_mw_wavelet_maps(L)
x = np.concatenate((x_scal_mw, x_wav_mw))

# perform inverse transform via hp
x_scal_hp_lm = lm2lm_hp(map2alm_mw(x_scal_mw, L, 0), L)
x_wav_hp_lm = np.zeros((L*(L+1)//2, nwvlts), dtype=np.complex)
offset = mw_size(L)
for i in range(nwvlts):
    lm = map2alm_mw(x_wav_mw[i*offset:(i+1)*offset], L, 0)
    x_wav_hp_lm[:,i] = lm2lm_hp(lm, L)
y_hp_lm = synthesis_axisym_lm_wav(x_wav_hp_lm, x_scal_hp_lm, B, L, J_min)
y = alm2map_mw(lm_hp2lm(y_hp_lm, L), L, 0).real

# adj: input = map, output = wavelet/scaling coeffs
# some random map
f = random_mw_map(L)

# perform inverse adjoint transform
g_wav_mw, g_scal_mw = synthesis_adjoint_axisym_wav_mw(f, B, L, J_min)
g = np.concatenate((g_scal_mw, g_wav_mw))

LHS = f.T.dot(y)
RHS = g.T.dot(x)
print(f"Synthesis: {LHS}, {RHS}")
print(f"Dot product test error: {LHS - RHS}")


##########
f1_lm = random_lms(L)
f2_lm = random_lms(L)

f1_lm_hp = lm2lm_hp(f1_lm, L)
f2_lm_hp = lm2lm_hp(f2_lm, L)

f1_wav_lm_hp, f1_scal_lm_hp = analysis_axisym_lm_wav(f1_lm_hp, B, L, J_min)
f2_wav_lm_hp, f2_scal_lm_hp = analysis_axisym_lm_wav(f2_lm_hp, B, L, J_min)

f1_scal = alm2map_mw(lm_hp2lm(f1_scal_lm_hp, L), L, 0)
f2_scal = alm2map_mw(lm_hp2lm(f2_scal_lm_hp, L), L, 0)
f1_wav = np.zeros(mw_size(L) * nwvlts, dtype=np.complex)
f2_wav = np.zeros(mw_size(L) * nwvlts, dtype=np.complex)
offset = mw_size(L)
for j in range(nwvlts):
    f1_wav[j*offset:(j+1)*offset] = alm2map_mw(
        lm_hp2lm(
            np.ascontiguousarray(f1_wav_lm_hp[:, j]), L
        ), L, 0
    )
    f2_wav[j*offset:(j+1)*offset] = alm2map_mw(
        lm_hp2lm(
            np.ascontiguousarray(f2_wav_lm_hp[:, j]), L
        ), L, 0
    )

f1 = alm2map_mw(f1_lm, L, 0)
f2_rec_adj = analysis_adjoint_axisym_wav_mw(f2_wav, f1_wav, B, L, J_min)

dot_error = f2_wav.dot(f1_wav) + f2_scal.dot(f1_scal) - f2_rec_adj.dot(f1)
print(f"Dot error: {dot_error}")