from pys2let import *
import numpy as np


def random_lms(L):
    return np.random.rand(L * L).astype(np.complex)


def random_mw_map(L, spin):
    return alm2map_mw(random_lms(L), L, spin)


def random_wavlet_maps(L, spin, nwvlts):
    f_scal = random_mw_map(L, spin)
    f_wav = np.concatenate([random_mw_map(L, spin) for _ in range(nwvlts)])
    return f_scal, f_wav


L = 10
B = 2
J_min = 2
N = 1
spin = 0
upsample = 1
nwvlts = pys2let_j_max(B, L, J_min) - J_min + 1


#### as on website http://sepwww.stanford.edu/sep/prof/pvi/conj/paper_html/node9.html
x = random_mw_map(L, spin)
y_scal, y_wav = random_wavlet_maps(L, spin, nwvlts)

y = analysis_adjoint_wav2px(y_wav, y_scal, B, L, J_min, N, spin, upsample)
x_wav, x_scal = analysis_px2wav(x, B, L, J_min, N, spin, upsample)

dot_product_error = y_wav.conj().dot(x_wav)
dot_product_error += y_scal.conj().dot(x_scal)
dot_product_error -= y.conj().dot(x)
print(f"ANALYSIS Dot product error: {dot_product_error}")

x_scal, x_wav = random_wavlet_maps(L, spin, nwvlts)
y = random_mw_map(L, spin)

x = synthesis_wav2px(x_wav, x_scal, B, L, J_min, N, spin, upsample)
y_wav, y_scal = synthesis_adjoint_px2wav(y, B, L, J_min, N, spin, upsample)

dot_product_error = y.conj().dot(x)
dot_product_error -= y_wav.conj().dot(x_wav)
dot_product_error -= y_scal.conj().dot(x_scal)
print(f"SYNTHESIS Dot product error: {dot_product_error}")
