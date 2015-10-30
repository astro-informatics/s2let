% s2let_demo11_plotLightProbeMaps
% -----------------------------------------------------------
% Plot the 'light probe' maps from FITS file
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 11th Aug 2015
% -----------------------------------------------------------
% ssht_plot_mollweide
%  'ColourBar'       = { false        [do not add colour bar (default)],
%                        true         [add colour bar] }
%  'Mode'            = { 0            Plot amplitude, or modulus is f complex (default),
%                        1            Plot real part,
%                        2            Plot imaginaty part,
%                        3            Plot modulus and arrows for real/img angle } 

clear; % all;
% close all ;

% ---------------
% Set the reality and upsample flag for signal analysis:
% ---------------
reality  = true;    % true for real data;
upsample = true ;  % true for full resolution plot

% ---------------
% Load data set (light probe map): 
% ---------------
%filename= '../../../data/LightProbe_ringmaps/heal_grace_ring.fits'; %Grace Cathedral, San Francisco
filename= '../../../data/LightProbe_ringmaps/heal_stpeters_ring.fits'; % St. Perter's Basilica, Rome
%filename= '../../../data/LightProbe_ringmaps/heal_galileo_ring.fits'; %The Galileo's Tomb, Florence
% filename= '../../../data/LightProbe_ringmaps/heal_rnl_ring.fits';  %Eucalyptus Grove, UC Berkeley
%filename= '../../../data/LightProbe_ringmaps/heal_uffizi_ring.fits';   %The Uffizi Gallery, Florence
% CMB: 
% filename = '../../../data/COM_CompMap_Synchrotron-commander_0256_R2.00.fits';
% filename = '../../../data/COM_CompMap_freefree-commander_0256_R2.00.fits';
% filename = '../../../data/COM_CompMap_CMB-commander_0256_R2.00.fits';
% filename = '../../../data/wmap_mcmc_base_k_synch_stk_u_9yr_v5.fits';
% filename =  '../../../data/somecmbsimu_hpx_128.fits'

% ----------
% Check what is in the FITS file
% ----------
% info =fitsinfo(filename);
% disp(info.Contents);
% info.BinaryTable
% rowend = info.BinaryTable.Rows;

% ----------
% Read and plot the FITS file
% ----------
[data_hpxmap, nside] = s2let_hpx_read_real_map(filename);

figure; 
s2let_hpx_plot_mollweide(data_hpxmap);
colorbar;
v = caxis;
temp = max(abs(v))*0.07; 
caxis([0  temp])
colormap bone
         
% ----------
% Convert from hpx map to mw map
% ----------
f_mw = s2let_hpx2mw(data_hpxmap);
% ----------
% Plot 
% ----------
mode = 1;
Lguessed= 2*nside;
% clip data 
max_f_mw = 1.     
f_mw(f_mw > max_f_mw) = max_f_mw;
figure; 
ssht_plot_mollweide(f_mw, Lguessed, 'Mode', mode);
colorbar;
v = caxis;
temp = max(abs(v))*0.1; 
caxis([0  temp])
colormap hot
colormap bone
% ----------
% Bandlimit the data
% ----------
L = 256; % 512   
B = 2;
J_min = 2;
nside_recon = L/2; %256 
flm = s2let_hpx_map2alm(data_hpxmap, 'L', L);
f = s2let_hpx_alm2map(flm, nside_recon, 'L', L);
% ----------
% Convert from hpx map to mw map
% ----------
f_bandlimit_mw = s2let_hpx2mw(f);

figure; 
ssht_plot_mollweide(f_bandlimit_mw, L, 'Mode', mode);
colorbar;
v = caxis;
temp = max(abs(v)); 
caxis([0  temp])
colormap bone

% ----------
% Clip data
% ----------
max_f_bandlimit_mw = 1.     
f_bandlimit_mw(f_bandlimit_mw > max_f_bandlimit_mw) = max_f_bandlimit_mw;
min_f_bandlimit_mw = -1.     
f_bandlimit_mw(f_bandlimit_mw < min_f_bandlimit_mw) = min_f_bandlimit_mw;
% ----------
% Plot 
% ----------
mode = 1;
figure; 
ssht_plot_mollweide(f_bandlimit_mw, L, 'Mode', mode);
colorbar;
v = caxis;
temp = max(abs(v)); 
caxis([0  temp])
colormap bone


