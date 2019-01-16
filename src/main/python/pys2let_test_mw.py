
import healpy as hp
from pys2let import *
import math
import matplotlib.pyplot as plt
import os

L = 32
spin = 0

# fname = '/Users/bl/Dropbox/Astrodata/SDSS/Fields/Planck_EBV_256rQ.fits'
fname = os.path.join(os.path.dirname(__file__), os.pardir,
                     os.pardir, os.pardir, 'data', 'somecmbsimu_hpx_128.fits')

# Read healpix map and compute alms.
# f_lm has size L*(L+1)/2
f_hpx = hp.read_map(fname)  # Initial map
f_lm = hp.map2alm(f_hpx, lmax=L - 1)  # Its alms

# Convert to MW sampling from spherical harmonics
f_mw = alm2map_mw(f_lm, L, spin)
f_lm_rec = map2alm_mw(f_mw, L, spin)
f_mw_rec = alm2map_mw(f_lm, L, spin)

Lpix = 3 * L

f_mw = alm2map_mw(f_lm, L, spin)
f_lm_rec = map2alm_mw(f_mw, L, spin)
f_mw_rec = alm2map_mw(f_lm, L, spin)


def zeropad(flm, Lin, Lout):
    f_lm_out = np.zeros([int(Lout * (Lout + 1) / 2), ], dtype=complex)
    for el in range(Lin):
        for em in range(el + 1):
            f_lm_out[healpy_lm(el, em, Lout)] = f_lm[healpy_lm(el, em, Lin)]
    return f_lm_out


f_lm_ext = zeropad(f_lm, L, Lpix)
f_mw_ext = alm2map_mw(f_lm_ext, Lpix, spin)
# f_mw_ext2 = alm2map_mw_bl(f_lm, 0, L, Lpix, spin)
# f_lm_ext2 = map2alm_mw_bl(f_mw_ext2, 0, L, Lpix, spin)
# f_mw_rec2 = alm2map_mw(f_lm_ext2, L, spin)

# Home made plotting routine! inputs : function f (1D array of MW signal), bandlimit L, plot axis ax, and title


def myplot(f, L, ax, title=''):
    # Compute the theta and phi values of the MW equiangular grid.
    thetas, phis = mw_sampling(L)
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi))
    ax.imshow(np.real(arr), vmin=0, vmax=1)
    if L > 10:
        step_phi = int(L / 4)
        step_theta = int(L / 4)
    else:
        step_phi = 3
        step_theta = 3
    selec = np.arange(0, nphi, step_phi)  # subset of phi tick labels
    ax.set_xticks(selec)
    ax.set_xticklabels(['%.1f' % x for x in phis[selec]])
    ax.set_xlabel(r'$\phi$')
    selec = np.arange(0, ntheta, step_theta)  # subset of theta tick labels
    ax.set_yticks(selec)
    ax.set_yticklabels(['%.1f' % x for x in thetas[selec]])
    ax.set_ylabel(r'$\theta$')
    ax.set_title(title)


# Plot equiangular map
fig, axs = plt.subplots(1, 2)
axs = axs.ravel()
myplot(f_mw, L, axs[0], 'L')
myplot(f_mw_ext, Lpix, axs[1], 'Lpix')
# myplot(f_mw_ext2, Lpix, axs[2], 'Lpix')
# myplot(f_mw_rec2, L, axs[3], 'L')
plt.show()
