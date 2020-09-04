from pys2let import *
import numpy as np

def random_lms(L):
    return np.random.rand(L * L).astype(np.complex)

def random_mw_map(L):
    return alm2map_mw(random_lms(L), L, 0)

def random_wavlet_maps(L, nwvlts):
    f_scal = random_mw_map(L)
    f_wav = np.concatenate([random_mw_map(L) for _ in range(nwvlts)])
    return f_scal, f_wav

L = 10
B = 2
J_min = 2
nwvlts = pys2let_j_max(B, L, J_min) - J_min + 1

#### as in s2let_test.c
f = random_mw_map(L)
f_rec = random_mw_map(L)

f_wav, f_scal = analysis_axisym_wav_mw(f, B, L, J_min)
f_wav_rec, f_scal_rec = analysis_axisym_wav_mw(f_rec, B, L, J_min)

f_rec = analysis_adjoint_axisym_wav_mw(f_wav, f_scal, B, L, J_min)

dot_product_error = f_wav_rec.conj().dot(f_wav)
dot_product_error += f_scal_rec.conj().dot(f_scal)
dot_product_error -= f_rec.conj().dot(f)

print(f"Dot product error: {dot_product_error}")

#### as on website http://sepwww.stanford.edu/sep/prof/pvi/conj/paper_html/node9.html
x = random_mw_map(L)
y_scal, y_wav = random_wavlet_maps(L, nwvlts)

y = analysis_adjoint_axisym_wav_mw(y_wav, y_scal, B, L, J_min)
x_wav, x_scal = analysis_axisym_wav_mw(x, B, L, J_min)

dot_product_error = y_wav.conj().dot(x_wav)
dot_product_error += y_scal.conj().dot(x_scal)
dot_product_error -= y.conj().dot(x)
print(f"ANALYSIS Dot product error: {dot_product_error}")

x_scal, x_wav = random_wavlet_maps(L, nwvlts)
y = random_mw_map(L)

x = synthesis_axisym_wav_mw(x_wav, x_scal, B, L, J_min)
y_wav, y_scal = synthesis_adjoint_axisym_wav_mw(y, B, L, J_min)

dot_product_error = y.conj().dot(x)
dot_product_error -= y_wav.conj().dot(x_wav)
dot_product_error -= y_scal.conj().dot(x_scal)
print(f"SYNTHESIS Dot product error: {dot_product_error}")
