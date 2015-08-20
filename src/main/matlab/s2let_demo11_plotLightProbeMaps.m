% s2let_demo11_plotLightProbeMaps
% -----------------------------------------------------------
% Plot the 'light probe' maps from FITS file
% -----------------------------------------------------------
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% 
% Modified S2LET package to perform Curvelets on the Sphere.
% -----------------------------------------------------------
% L is the angular band-limit.
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 11th Aug 2015
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
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
reality = false;    % true for real data;
upsample = true ;  % true for full resolution plot

% ---------------
% Load data set (light probe map): 
% ---------------
%filename= '../../../data/heal_grace.fits'; %Grace Cathedral, San Francisco (dynamical range 200, 0000:1) 
%filename= '../../../data/heal_stpeters.fits'; % St. Perter's Basilica, Rome (dynamical range 200, 0000:1) 
%filename= '../../../data/heal_galileo.fits'; %The Galileo's Tomb, Florence (dynamical range 7000:1 original)
%filename= '../../../data/heal_rnl.fits';  %Eucalyptus Grove, UC Berkeley (dynamical range 5000:1 original)
filename= '../../../data/heal_uffizi.fits';   %The Uffizi Gallery, Florence (dynamical range 500:1 original)
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

% datacell = fitsread(filename,'binarytable');
% data = datacell{1};
% sz = size(data);

% Complex:

f_hpx = fitsread(filename,'binarytable',...
               'TableColumns',[1],...
               'TableRows', [256])     

flm_gen= cell2mat(f_hpx);
L = 32;
f_gen = ssht_inverse(flm_gen, L,'Reality', false);
ssht_plot_mollweide(f_gen, L);
colorbar; 
%}

% ----------
% Read and plot the FITS file
% ----------
[data_hpxmap, nside] = s2let_hpx_read_real_map(filename);
% figure; 
% s2let_hpx_plot_mollweide(data_hpxmap);
sz = size(data_hpxmap);
Lguessed = 2*nside 
% whos hpxmap                                                     
                                                            

% ----------
% Convert from hpx map to mw map
% ----------
f_mw = s2let_hpx2mw(data_hpxmap);
% whos data_hpxmap
% whos f_mw
% ----------
% Plot 
% ----------
mode = 1;
figure; 
ssht_plot_mollweide(f_mw, Lguessed, 'Mode', mode);
colorbar;
v = caxis;
temp = max(abs(v)); 
caxis([0  temp])

% ----------
% bandlimit the data
% ----------
L = 128; % 512
B = 2;
J_min = 2;
nside_recon = 64; %256 
flm = s2let_hpx_map2alm(data_hpxmap, 'L', L);
f = s2let_hpx_alm2map(flm, nside_recon, 'L', L);
% ----------
% Convert from hpx map to mw map
% ----------
f_bandlimit_mw = s2let_hpx2mw(f);
% ----------
% Plot 
% ----------
mode = 1;
figure; 
ssht_plot_mollweide(f_bandlimit_mw, L, 'Mode', mode);
colorbar;
v = caxis;
temp = max(abs(v));
caxis([0 temp])


