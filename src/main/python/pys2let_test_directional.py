
import healpy as hp
from pys2let import *
import math
import matplotlib.pyplot as plt

nside = 512
L = 16
J_min = 4
B = 1.5 #DOUBLE
N = 3 # Number of directions
spin = 0 # are we dealing with spin signals? set to 0 for temperature. if non-zero, plotting routines must be changed!
upsample = 0 # 1 means all scales at full resolution L # 0 means multiresolution wavelet transform

scal_lmax = B**J_min
J = pys2let_j_max(B, L, J_min) # Compute maximum scale
print 'Jmax = ', J

# The filename of some random healpix map
fname = '/Users/keir/Documents/s2let_ilc/simu_data/camb_simu_map.fits'
#fname = '/Users/bl/Dropbox/Astrodata/SDSS/Fields/Planck_EBV_256rQ.fits'
#fname = '/Users/bl/Dropbox/Wavelets/s2let/data/somecmbsimu_hpx_128.fits'

# Read healpix map and compute alms.
# f_lm has size L*(L+1)/2
f_ini = hp.read_map(fname) # Initial map
f_lm = hp.map2alm(f_ini, lmax=L-1) # Its alms
f = hp.alm2map(f_lm, nside=nside, lmax=L-1) # Band limited version

# Call pys2let and compute directional wavelet transform.
# Returns an array of MW maps
# The way to access them is described below.
print 'Running analysis_lm2wav'
f_wav, f_scal = analysis_lm2wav(f_lm, B, L, J_min, N, spin, upsample)
print 'Done'


print 'size f_scal f_wav', f_scal.size, f_wav.size, f_wav.size / f_scal.size

# Uses synthesis to reconstruct the input alms.
print 'Running synthesis_wav2lm'
f_lm_rec = synthesis_wav2lm(f_wav, f_scal, B, L, J_min, N, spin, upsample)
print 'Done'
f_rec = hp.alm2map(f_lm_rec, nside=nside, lmax=L-1)
#plot to compare the input/output maps
#hp.mollview(f)
#hp.mollview(f_rec)

# Convert to MW sampling from spherical harmonics
f_mw = alm2map_mw(f_lm, L, spin)
f_lm_rec = map2alm_mw(f_mw, L, spin) # To check that the analysis routine works
f_mw_rec = alm2map_mw(f_lm, L, spin)


# Home made plotting routine! inputs : function f (1D array of MW signal), bandlimit L, plot axis ax, and title
def myplot(f, L, ax, title=''):
    thetas, phis = mw_sampling(L) # Compute the theta and phi values of the MW equiangular grid.
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi)) # Convert the input MW 1D array into 2D array with rows for theta and columns for phi. As simple as that!
    ax.imshow(arr.astype(float),interpolation='none') # Don't forget to convert to float, otherwise imshow will complain about complex numbers. For polalisation we will have to look into other operations (real part, imag part, angle).

    # This is to undersample the grid of x/yticklabels.
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
'''fig, ax = plt.subplots(1,1)
myplot(f_mw, L, ax, 'Input map converted to MW')
fig, ax = plt.subplots(1,1)
myplot(f_mw_rec, L, ax, 'Input map converted to MW (reconstructed)')
fig.savefig('test_directional_python_wrappers_1.png')'''

# Plot scaling function map
fig, ax = plt.subplots(1,1)
myplot(f_scal, scal_lmax, ax, title='Scaling function map') #NEEDS CORRECTING

# Create giant array figure
fig, axs = plt.subplots(J-J_min+1, N, figsize=(4*N, 3*(J-J_min)))
axs = axs.ravel()
# Loop through scales j and directions n
for j in range(J_min, J+1):
    for n in range(0, N):
        # Retreive the boundaries and positions of the right wavelet scale in the giant f_wav array!
        offset, bandlimit, nelem, nelem_wav = wav_ind(j, n, B, L, N, J_min, upsample)
        # The right wavelet map corresponding to (j,n) will be f_wav[offset:offset+nelem].
        # It has a band-limit bandlimit
        # nelem_wav is the total number of elements in the j-th scale (i.e., sum of all directions). It's a safety check to verify that we are not forgetting any directions.
        print 'plot id', (j-J_min)*N+n, 'j=', j, 'n=', n, 'bounds=',offset, 'to', offset+nelem, 'Total elems:', nelem_wav
        # Make the plot!
        myplot(f_wav[offset:offset+nelem], bandlimit, axs[(j-J_min)*N+n],
            title='Scale '+str(j+1-J_min)+'/'+str(J-J_min+1)+', direction '+str(n+1)+'/'+str(N))

# Pretty adjustment
fig.subplots_adjust(hspace=0.4, wspace=0.5)
fig.savefig('test_directional_python_wrappers_2.png')
plt.show()
