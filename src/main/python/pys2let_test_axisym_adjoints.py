from pys2let import *
import numpy as np

# DOT PRODUCT TESTS
# if y = Ax and g = A'f, show that f'y = g'x
L = 10
B = 2
J_min = 2
J = pys2let_j_max(B, L, J_min)
nwvlts = J - J_min + 1

def random_mw_map(complex=False):
    f = np.random.rand(mw_size(L))
    if complex:
        f = f + np.random.rand(mw_size(L)) * 1j
    return f

def random_mw_wavelet_maps(complex=False):
    f_wav_mw = np.random.rand(mw_size(L) * nwvlts)
    f_scal_mw = np.random.rand(mw_size(L))
    if complex:
        f_wav_mw = f_wav_mw + np.random.rand(mw_size(L) * nwvlts) * 1j
        f_scal_mw = f_scal_mw + np.random.rand(mw_size(L)) * 1j
    return f_scal_mw, f_wav_mw


#  ANALYSIS
complex = False
# fwd: input = map, output = wavelet/scaling coeffs
# some random input map
x = random_mw_map(complex)  # MW map

# perform fwd transform via hp
x_lm = map2alm_mw(x.astype(np.complex), L, 0)
x_hp_lm = lm2lm_hp(x_lm, L)
y_wav_hp_lm, y_scal_hp_lm = analysis_axisym_lm_wav(x_hp_lm, B, L, J_min)
y_scal_mw = alm2map_mw(lm_hp2lm(y_scal_hp_lm, L), L, 0)
y_wav_mw = np.zeros((mw_size(L), nwvlts), dtype=np.complex)
for i in range(nwvlts):
    y_wav_mw[:,i] = alm2map_mw(lm_hp2lm(np.ascontiguousarray(y_wav_hp_lm[:,i]), L), L, 0)
y = np.concatenate((y_scal_mw, y_wav_mw.flatten("F"))).real  # MW 

# adj: input = wavelet/scaling coeffs, output = map
# some random wavelet coeffs
f_scal_mw, f_wav_mw = random_mw_wavelet_maps(complex)  
f = np.concatenate((f_scal_mw, f_wav_mw))  # MW 

# perform adjoint transform
g = analysis_adjoint_axisym_wav_mw(f_wav_mw, f_scal_mw, B, L, J_min)  # MW map

LHS = f.T.dot(y)
RHS = g.T.dot(x)
print(f"Analysis: {LHS}, {RHS}")
print(f"Dot product test error: {LHS - RHS}")


# SYNTHESIS
complex = False
# inv: input = wavelet/scaling coeffs, output = map 
# some random wavelet coeffs
x_scal_mw, x_wav_mw = random_mw_wavelet_maps(complex)
x = np.concatenate((x_scal_mw, x_wav_mw))

# perform inverse transform via hp
x_scal_hp_lm = lm2lm_hp(map2alm_mw(x_scal_mw.astype(np.complex), L, 0), L)
x_wav_hp_lm = np.zeros((L*(L+1)//2, nwvlts), dtype=np.complex)
offset = mw_size(L)
for i in range(nwvlts):
    lm = map2alm_mw(x_wav_mw[i*offset:(i+1)*offset].astype(np.complex), L, 0)
    x_wav_hp_lm[:,i] = lm2lm_hp(lm, L)
y_hp_lm = synthesis_axisym_lm_wav(x_wav_hp_lm, x_scal_hp_lm, B, L, J_min)
y = alm2map_mw(lm_hp2lm(y_hp_lm, L), L, 0).real

# adj: input = map, output = wavelet/scaling coeffs
# some random map
f = random_mw_map(complex)

# perform inverse adjoint transform
g_wav_mw, g_scal_mw = synthesis_adjoint_axisym_wav_mw(f, B, L, J_min)
g = np.concatenate((g_scal_mw, g_wav_mw))

LHS = f.T.dot(y)
RHS = g.T.dot(x)
print(f"Synthesis: {LHS}, {RHS}")
print(f"Dot product test error: {LHS - RHS}")