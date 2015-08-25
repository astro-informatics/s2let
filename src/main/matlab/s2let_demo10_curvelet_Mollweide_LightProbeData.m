% s2let_demo10_curvelet_Mollweide_LightProbeData 
%
% Plot curvelet coefficients on multiple Mollweide projections.
% The function generates one plot of the scaling function
% contribution and a grid of plots for each orientation of
% each scale of the curvelet contributions. 
% -----------------------------------------------------------
% f_cur is cell array with all the curvelet coefficients.
% Its first index is the curvelet scale j, the second
% index is the orientation g, and each element is a
% function on the sphere in MW sampling.
% scal is the corresponding scaling function contribution
% (i.e. just a single function on the sphere).
%
% B is the wavelet parameter.
% L is the angular band-limit.
% N is the orientational band-limit.
% J_min is the first curvelet scale in f_cur.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%  'Spin'            = { non-negative integers (default=0) }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'Upsample'        = { false        [multiresolution algorithm (default)],
%                        true         [full resolution curvelets] }
%  'SpinLowered'     = { true         [Apply normalisation factors for 
%                                      spin-lowered curvelets and 
%                                      scaling function.],
%                        false        [Apply the usual normalisation factors such
%                                      that the curvelets fulfill the admissibility
%                                      condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the curvelets
%                       should be lowered from (default = 0)]
%
%
% For ssht_plot_mollweide(f, L, <options>), 
% f is the sampled function and L is the harmonic band-limit.
% Valid options include: 
%  'ColourBar'       = { false        [do not add colour bar (default)],
%                        true         [add colour bar] }
%  'Mode'            = { 0            Plot amplitude, or modulus is f complex (default),
%                        1            Plot real part,
%                        2            Plot imaginaty part,
%                        3            Plot modulus and arrows for real/img angle }
%  'Spin'            = { non-negative integers (default=0) }
% -----------------------------------------------------------
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt, Martin Büttner,
%                     Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------
 
clear; % all;
% close all ;

% ---------------
% Load data set (light probe map): 
% ---------------
%filename= '../../../data/LightProbe_ringmaps/heal_grace_ring.fits';    % Grace Cathedral, San Francisco  
%filename= '../../../data/LightProbe_ringmaps/heal_stpeters_ring.fits'; % St. Perter's Basilica, Rome 
%filename= '../../../data/LightProbe_ringmaps/heal_galileo_ring.fits';  % The Galileo's Tomb, Florence
% filename= '../../../data/LightProbe_ringmaps/heal_rnl_ring.fits';     % Eucalyptus Grove, UC Berkeley 
filename= '../../../data/LightProbe_ringmaps/heal_uffizi_ring.fits';    % The Uffizi Gallery, Florence  
% ----------
% Read data from the FITS file
% ----------
[data_hpxmap, nside] = s2let_hpx_read_real_map(filename);

% ----------
% Band-limit the data
% ----------
L = 256;             % 512
nside_recon = L/2;  % 256 
flm = s2let_hpx_map2alm(data_hpxmap, 'L', L);
f = s2let_hpx_alm2map(flm, nside_recon, 'L', L);

% ----------
% Convert from hpx map to mw map
% ----------
f_bandlimit_mw = s2let_hpx2mw(f);

% ----------
% Clip data
% ----------
max_f_bandlimit_mw = 1.     
f_bandlimit_mw(f_bandlimit_mw > max_f_bandlimit_mw) = max_f_bandlimit_mw;

% ---------------
% Define curvelet parameters: 
% ---------------
B = 2;   
Spin = 0; 
N= L;      
J_min = 3;  
J =s2let_jmax(L, B);  
 
% ---------------
% Define plotting parameters:
% ---------------
zoomfactor =1.;  
ny = 16; 
nx = 4;  

maxfigs = nx*ny;
pltroot = '../../../figs' ; 
configstr = ['N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]; 

% ---------------
% Set the reality and upsample flag for signal analysis : 
% ---------------
reality = true;    % true for real data;
upsample = true ;  % true for full resolution plot 

% -----------------
% Signal analysis: 
% -----------------
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_bandlimit_mw,  ...
                                                           'B', B, 'L', L, ...
                                                           'J_min', J_min, ...
                                                           'Spin', Spin, ...
                                                           'Reality', reality, ...
                                                           'Upsample', upsample , ...
                                                           'SpinLowered', false, ...
                                                           'SpinLoweredFrom', 0);


figure('Position',[20 20 1700 1400]) %100 100 1300 1000
subplot(ny, nx, 1);
% --- plot initial data --- % 
ssht_plot_mollweide(f_bandlimit_mw, L);
%
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])

%
subplot(ny, nx, 2);
% --- plot scaling function contributions --- %
if (upsample ~= false) %full resolution 
    ssht_plot_mollweide(f_scal, L);
else
    band_limit = min([ s2let_bandlimit(J_min-1,J_min,B,L) L ]);
    ssht_plot_mollweide(f_scal, band_limit);
end
%
title('Scaling fct')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])


% --- plot curvelet kernel contributions --- % 
ind = 2;
for j = J_min:J  
    band_limit = min([ s2let_bandlimit(j,J_min,B,L) L ]);  
    Nj = band_limit; 
    if (reality == false) % default i.e. complex signals
        enmax = 2*Nj-1; 
    else % i.e. real signals
        enmax = Nj; 
    end 
	for en = 1: enmax
		ind = ind + 1;
        if ind <= maxfigs
            subplot(ny, nx, ind);
            %
            if (upsample ~= false)
                ssht_plot_mollweide(reshape(f_cur{j-J_min+1}(en,:,:), L, 2*L-1), L);
            else
                ssht_plot_mollweide(reshape(f_cur{j-J_min+1}(en,:), band_limit, 2*band_limit-1), band_limit);
            end
            %
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Curvelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)],'FontSize', 10)
        end
    end
end
colormap(jet)
if (upsample ~= false)
    fname = [pltroot,'/s2let_demo10_', configstr, '_spin_curvelet_LightProbeMap_FullRes.png']
else
    fname = [pltroot,'/s2let_demo10_', configstr, '_spin_curvelet_LightProbeMap_MultiRes.png']
end
print('-r200', '-dpng', fname)

% ---------- 
% Compare the reconstructed signal with the initial signals: 
% ---------- 
%{
f_rec = s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                  'B', B, 'L', L, ...
                                                  'J_min', J_min, ...
                                                  'Spin', Spin, ...
                                                  'Reality', reality, ...
                                                  'Upsample', upsample, ...
                                                  'SpinLowered', false, ...
                                                  'SpinLoweredFrom', 0);

figure('Position',[100 100 900 200]) 
subplot(2, 2, 1);
ssht_plot_mollweide(f_gen, L);
title('initial signal')
hold on
subplot(2, 2, 2);
ssht_plot_mollweide(f_rec,L);
title('reconstructed signal')
if (upsample == false)
    fname = [pltroot,'/s2let_demo10_', configstr, '_spin_curvelet_EarthTomo_FullRes_Int_Rec_signal.png']
else
    fname = [pltroot,'/s2let_demo10_', configstr, '_spin_curvelet_EarthTomo_MultiRes_Int_Rec_signal.png']
end
print('-r200', '-dpng', fname)
% Check error:
check_error = max(abs(f_gen(:)-f_rec(:)))
%}                           

% ============
% Using directional wavelets for analysis 
% ============
%{
[f_wav, f_scal] = s2let_transform_analysis_mw(f_gen, 'B', B, 'J_min', J_min, ...
                                              'N', N, 'Reality',reality,...
                                              'Upsample', upsample, 'Spin', 0);

% FULL-RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(ny, nx, 1);
ssht_plot_mollweide(f_gen, L);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(ny, nx, 2);
ssht_plot_mollweide(f_scal, L);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
title('Scaling fct')
ind = 2;
for j = J_min:J
	for en = 1:N
		ind = ind + 1;
        if ind <= maxfigs
            subplot(ny, nx, ind);
            ssht_plot_mollweide(f_wav{j-J_min+1,en}, L, 'Mode', 1);
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Wavelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)])
        end
	end
end
colormap(jet)
fname = [pltroot,'/s2let_demo10_', configstr, '_directional_wavelets_LightProbeMap_FullRes.png']
print('-r200', '-dpng', fname)
%}