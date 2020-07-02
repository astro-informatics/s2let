from pys2let import *
import numpy as np
import healpy as hp

# DOT PRODUCT TESTS
# if y = Ax and g = A'f, show that f'y = g'x
Nside = 32
L = 10
B = 2
J_min = 2
J = pys2let_j_max(B, L, J_min)
nwvlts = J - J_min + 1

def random_mw_map():
    return np.random.rand(mw_size(L))

def random_mw_wavelet_maps():
    f_wav_mw = np.random.rand(mw_size(L) * nwvlts)
    f_scal_mw = np.random.rand(mw_size(L))
    return f_scal_mw, f_wav_mw


#  ANALYSIS
# fwd: input = map, output = wavelet/scaling coeffs
# some random input map
x = random_mw_map()  # MW map

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
f_scal_mw, f_wav_mw = random_mw_wavelet_maps()  
f = np.concatenate((f_scal_mw, f_wav_mw))  # MW 

# perform adjoint transform
g = analysis_adjoint_axisym_wav_mw(f_wav_mw, f_scal_mw, B, L, J_min)  # MW map

LHS = f.T.dot(y)
RHS = g.T.dot(x)
print(f"Analysis: {LHS}, {RHS}")
