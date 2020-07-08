from pys2let import *
import numpy as np

def random_lms(L):
    return np.random.rand(L * L).astype(np.complex)  # not right for N > 1

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

#### as in s2let_test.c
f = random_mw_map(L, spin)
f_rec = random_mw_map(L, spin)

f_wav, f_scal = analysis_px2wav(f, B, L, J_min, N, spin, upsample)
f_wav_rec, f_scal_rec = analysis_px2wav(f_rec, B, L, J_min, N, spin, upsample)

f_rec = analysis_adjoint_wav2px(f_wav, f_scal, B, L, J_min, N, spin, upsample)

dot_product_error = f_wav_rec.conj().dot(f_wav)
dot_product_error += f_scal_rec.conj().dot(f_scal)
dot_product_error -= f_rec.conj().dot(f)

print(f"Dot product error: {dot_product_error}")

#### as on website http://sepwww.stanford.edu/sep/prof/pvi/conj/paper_html/node9.html
x = random_mw_map(L, spin)
y_scal, y_wav = random_wavlet_maps(L, spin, nwvlts)

y = analysis_adjoint_wav2px(y_wav, y_scal, B, L, J_min, N, spin, upsample)
x_wav, x_scal = analysis_px2wav(x, B, L, J_min, N, spin, upsample)

dot_product_error = y_wav.conj().dot(x_wav)
dot_product_error += y_scal.conj().dot(x_scal)
dot_product_error -= y.conj().dot(x)
print(f"Dot product error: {dot_product_error}")
