import healpy as hp
from pys2let import *
import math
import matplotlib.pyplot as plt
import os

nside = 128
L = 128
J_min = 1
B = 3
J = pys2let_j_max(B, L, J_min)

# The filename of some random healpix map
fname = os.path.join(os.path.dirname(__file__), os.pardir,
                     os.pardir, os.pardir, 'data', 'somecmbsimu_hpx_128.fits')

# Read healpix map and compute alms.
# f_lm has size L*(L+1)/2
f_ini = hp.read_map(fname)  # Initial map
f_lm = hp.map2alm(f_ini, lmax=L - 1)  # Its alms
f = hp.alm2map(f_lm, nside=nside, lmax=L - 1)  # Band limited version

hp.mollview(f)

# Call pys2let and compute wavelet transform. Returns the harmonic coefficients of the wavelets.
# f_scal_lm has size L*(L+1)/2
# f_wav_lm has size L*(L+1)/2 by J-J_min+1
f_wav_lm, f_scal_lm = analysis_axisym_lm_wav(f_lm, B, L, J_min)

# Reconstruct healpix maps on the sphere and plot them
f_scal = hp.alm2map(f_scal_lm, nside=nside, lmax=L - 1)
hp.mollview(f_scal)
f_wav = np.empty([12 * nside * nside, J - J_min + 1])
for j in range(J - J_min + 1):
    flm = f_wav_lm[:, j].ravel()
    f_wav[:, j] = hp.alm2map(flm, nside=nside, lmax=L - 1)
    hp.mollview(f_wav[:, j])

# Uses synthesis to reconstruct the input map.
f_lm_rec = synthesis_axisym_lm_wav(f_wav_lm, f_scal_lm, B, L, J_min)
f_rec = hp.alm2map(f_lm_rec, nside=nside, lmax=L - 1)
hp.mollview(f_rec)

plt.show()
