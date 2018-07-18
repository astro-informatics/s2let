import healpy as hp
from pys2let import *
import math
import matplotlib.pyplot as plt
from copy import copy



L = 3024
J_min = 3

Bs = np.array([2, 1.4, 1.2, 1.1])
L_transitions = np.array([200, 1100, 1900])
hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits, J, L_bounds = construct_hybrid_tiling(L, J_min, Bs, L_transitions)

print 'Bs', Bs
print 'L_bounds', L_bounds
fig, axs = plt.subplots(Bs.size+2, 1, figsize=(8, 14))
axs = axs.ravel()
J_mins = np.zeros((Bs.size,), dtype=np.int32)
J_mins[0] = J_min
for k in range(Bs.size):
    scal_l, wav_l = axisym_wav_l(Bs[k], L, J_mins[k])
    Jt = pys2let_j_max(Bs[k], L, 0)
    axs[k].plot(scal_l, lw=2, ls='dashed')
    for j in range(0,Jt+1-J_mins[k]):
        axs[k].plot(wav_l[:,j], lw=2, ls='solid')
    for kk in range(Bs.size):
        axs[k].axvline(L_bounds[kk], c='k', ls='dashed')
    axs[k].set_xscale('log')
    axs[k].set_ylim([0, 1.2])
    axs[k].set_xlim([1, L])
    axs[k].set_title('Tiling %i'%(k+1)+' on %i'%Bs.size+', with B=%.2f'%Bs[k]+', Jmin=0, and defined at %i'%L_bounds[k]+'<el<%i'%L_bounds[k+1])
for k in [Bs.size, Bs.size+1]:
    axs[k].plot(hybrid_scal_l, lw=2, ls='dashed')
    for j in range(0,J+1):
        axs[k].plot(hybrid_wav_l[:,j], lw=2, ls='solid')
    for kk in range(Bs.size):
        axs[k].axvline(L_bounds[kk], c='k', ls='dashed')
    axs[k].set_ylim([0, 1.2])
    axs[k].set_xlim([1, L])
    if k == Bs.size:
        axs[k].set_title('Hybrid tiling (log scale)')
        axs[k].set_xscale('log')
    else:
        axs[k].set_title('Hybrid tiling (linear scale)')
plt.tight_layout()
#plt.show()

hybrid_wav_l = hybrid_wav_l.T.ravel()
res = verify_tiling(L, hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits)
print 'Is the tiling OK and giving an invertible wavelet transform?', res
fig.savefig('Hybrid_tiling.png')
#plt.show()
#stop


J_min = 2
nside = 64
N = 3
spin = 0
Bs = np.array([4, 2])
J_min = 0
L_transitions = np.array([64])
L = 256
hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits, J, L_bounds = construct_hybrid_tiling(L, J_min, Bs, L_transitions)
hybrid_wav_l = hybrid_wav_l.T.ravel()
res = verify_tiling(L, hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits)
print 'Is the tiling OK and giving an invertible wavelet transform?', res


# The filename of some arbitrary healpix map
fname = '/Users/bl/Dropbox/Astrodata/SDSS/Fields/Planck_EBV_256rQ.fits'
#fname = '/Users/bl/Dropbox/Wavelets/s2let/data/somecmbsimu_hpx_128.fits'

f_ini = hp.read_map(fname) # Initial map
f_lm = hp.map2alm(f_ini, lmax=L-1) # Its alms
f = hp.alm2map(f_lm, nside=nside, lmax=L-1) # Band limited version

print 'Running analysis_lm2wav'
f_wav, f_scal = analysis_lm2wav_manualtiling(f_lm, L, N, spin, hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits)

# Uses synthesis to reconstruct the input alms.
print 'Running synthesis_wav2lm'
f_lm_rec = synthesis_wav2lm_manualtiling(f_wav, f_scal, L, N, spin, hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits)
f_rec = hp.alm2map(f_lm_rec, nside=nside, lmax=L-1)

# Convert to MW sampling from spherical harmonics
f_mw = alm2map_mw(lm_hp2lm(f_lm, L), L, spin)
#f_lm_rec = map2alm_mw(f_mw, L, spin) # To check that the analysis routine works
#f_mw_rec = alm2map_mw(f_lm_rec, L, spin)
f_mw_rec = alm2map_mw(lm_hp2lm(f_lm_rec, L), L, spin)


print 'ACCURACY HARM:', np.max(f_mw_rec-f_mw), np.mean(f_mw_rec/f_mw), np.std(f_mw_rec/f_mw)
print 'ACCURACY REAL:', np.max(f_lm_rec-f_lm), np.mean(f_lm_rec/f_lm), np.std(f_lm_rec/f_lm)

# Home made plotting routine! inputs : function f (1D array of MW signal), bandlimit L, plot axis ax, and title
def myplot(f, L, ax, title=''):
    thetas, phis = mw_sampling(L)
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi))
    ax.imshow(arr.astype(float))
    if L > 10:
        step_phi = int(L/4)
        step_theta = int(L/4)
    else:
        step_phi = 3
        step_theta = 3
    selec = np.arange(0, nphi, step_phi) # subset of phi tick labels
    ax.set_xticks(selec)
    ax.set_xticklabels(['%.1f' % x for x in phis[selec]])
    ax.set_xlabel(r'$\phi$')
    selec = np.arange(0, ntheta, step_theta) # subset of theta tick labels
    ax.set_yticks(selec)
    ax.set_yticklabels(['%.1f' % x for x in thetas[selec]])
    ax.set_ylabel(r'$\theta$')
    ax.set_title(title)


# Plot equiangular map
fig, ax = plt.subplots(1,1)
myplot(f_mw, L, ax, 'Input map converted to MW')
fig, ax = plt.subplots(1,1)
myplot(f_mw_rec, L, ax, 'Input map converted to MW (reconstructed)')
fig.savefig('test_directional_python_wrappers_1.png')


# Create giant array figure
fig, axs = plt.subplots(J+1, N, figsize=(4*N, 3*(J)))
axs = axs.ravel()
# Loop through scales j and directions n
offset = 0
for j in range(0, J+1):
    for n in range(0, N):
        bandlimit = hybrid_wav_bandlimits[j]
        nelem = bandlimit*(2*bandlimit-1)
        #offset, bandlimit, nelem, nelem_wav = wav_ind(j, n, B, L, N, J_min, upsample)
        print 'plot id', (j)*N+n, 'j=', j, 'n=', n, 'bounds=',offset, 'to', offset+nelem
        # Make the plot!
        myplot(f_wav[offset:offset+nelem], bandlimit, axs[(j)*N+n],
            title='Scale '+str(j+1)+'/'+str(J+1)+', direction '+str(n+1)+'/'+str(N))
        offset += nelem

# Pretty adjustment
fig.subplots_adjust(hspace=0.4, wspace=0.5)
fig.savefig('test_directional_python_wrappers_2.png')
plt.show()
