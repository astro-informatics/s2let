% s2let_demo8 - curvelet analysis on the Earth topographic map
%
% Plot curvelet coefficients on multiple Mollweide projections.
% The function generates one plot of the scaling function
% contribution and a grid of plots for each orientation of
% each scale of the curvelet contributions.
%
% -----------------------------------------------------------
% f_cur is cell array with all the curvelet coefficients.
% Its first index is the curvelet scale j, the second
% index is the orientation g, and each element is a
% function on the sphere in MW sampling.
% scal is the corresponding scaling function contribution
% (i.e. just a single function on the sphere).
%
% B is the curvelet parameter. 
% L is the angular band-limit.
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
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

clear all;
close all ;

load('EGM2008_Topography_flms_L0128');
% ---------------
% band limit the signals: 
% --------------- 
L = 64;   % (set L>=64 to see multi-resolution effects) 
flm_gen = flm(1:L^2,1);
f_gen = ssht_inverse(flm_gen, L,'Reality', true);

% ---------------
% Set reality flag for the transform (to improve performance if f is real): 
% ---------------
reality = true;

% ---------------
% Define curvelet parameters: 
% ---------------
Spin = 0;             % Spin value of wavelet
B = 5;                % B=2 for dyadic sampling
J_min = 1;            % Minimum scale probed by curvelets 
J =s2let_jmax(L, B);  % Maximum scale probed by curvelets =ceil(log L/ log B); 
 
% ---------------
% Define Plotting parameters: 
% ---------------
zoomfactor= 1.0;  
ny = 7;  
nx = 3;  

maxfigs = nx*ny;
pltroot = '../../../figs' ; 
configstr = ['_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)];

% ================================================
% FULL RESOLUTION PLOT (Upsample: true)
% ================================================
% -----------------
% Signal analysis: 
% -----------------
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_gen,  ...
                                                           'B', B, 'L', L, ...
                                                           'J_min', J_min, ...
                                                           'Spin', Spin, ...
                                                           'Reality', reality, ...
                                                           'Upsample', true , ...
                                                           'SpinLowered', false, ...
                                                           'SpinLoweredFrom', 0);

figure('Position',[20 20 1500 1400]) %100 100 1300 1000
subplot(ny, nx, 1);
% --- plot initial data --- % 
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
%
title('Initial data','FontSize',8)
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
% 
subplot(ny, nx, 2);
% --- plot scaling function contributions --- % 
ssht_plot_mollweide(f_scal, L, 'Mode', 1);
%
title('Scaling fct','FontSize',8)
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
% --- plot curvelet kernel contributions --- % 
ind = 2;
for j = J_min:J  
    band_limit = min([ s2let_bandlimit(j,J_min,B,L) L ]);
    % Nj is the orientational band-limit at j-th scale.
    Nj = band_limit; 
    if (reality == false) 
        enmax = 2*Nj-1; 
    else                 
        enmax = Nj; 
    end 
	for en = 1: enmax
		ind = ind + 1;
        if ind <= maxfigs
            subplot(ny, nx, ind);
            %
            ssht_plot_mollweide(reshape(f_cur{j-J_min+1}(en,:,:), L, 2*L-1), L, 'Mode', 1);
            %
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Curvelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)],'FontSize', 8)
        end
    end
end
colormap(jet)
fname = [pltroot,'/s2let_demo8_', configstr, '_curvelet_EarthTopo_FullRes.png']
print('-r200', '-dpng', fname)


% ---------- 
% Plot reconstructed map (to compare with the initial map): 
% ---------- 
% -----------------
% Signal synthesis: 
% -----------------
f_rec = s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                  'B', B, 'L', L, ...
                                                  'J_min', J_min, ...
                                                  'Spin', Spin, ...
                                                  'Reality', reality, ...
                                                  'Upsample', true, ...
                                                  'SpinLowered', false, ...
                                                  'SpinLoweredFrom', 0);

figure('Position',[100 100 900 200]) 
subplot(2, 2, 1);
% --- plot initial data --- % 
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
title('initial signal')
hold on
subplot(2, 2, 2);
% --- plot reconstructed data --- % 
ssht_plot_mollweide(f_rec,L, 'Mode', 1);
title('reconstructed signal')
fname = [pltroot,'/s2let_demo8_', configstr, '_spin_curvelet_EarthTomo_FullRes_Int_Rec_signal.png']
print('-r200', '-dpng', fname)
% Check error:
check_error = max(abs(f_gen(:)-f_rec(:)))                   



% ================================================
% MULTI-RESOLUTION PLOT (Upsample: false)
% ================================================
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_gen,  ...
                                                           'B', B, 'L', L, ...
                                                           'J_min', J_min, ...
                                                           'Spin', Spin, ...
                                                           'Reality', reality, ...
                                                           'Upsample', false, ...
                                                           'SpinLowered', false, ...
                                                           'SpinLoweredFrom', 0);  

figure('Position',[20 20 1500 1400]) %100 100 1300 1000
subplot(ny, nx, 1);
% --- plot initial data --- % 
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
%
title('Initial data','FontSize',8)
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
% 
subplot(ny, nx, 2);
% --- plot scaling function contributions --- % 
bl = min([ s2let_bandlimit(J_min-1,J_min,B,L) L ]);
ssht_plot_mollweide(f_scal, bl, 'Mode', 1);
%
title('Scaling fct','FontSize',8)
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
% --- plot curvelet kernel contributions --- % 
ind = 2;
for j = J_min:J
    band_limit = min([ s2let_bandlimit(j,J_min,B,L) L ]);  
    Nj = band_limit; 
    if (reality == false) 
        enmax = 2*Nj-1; 
    else                 
        enmax = Nj; 
    end 
	for en = 1: enmax   
		ind = ind + 1;
        if ind <= maxfigs
            subplot(ny, nx, ind);
            %
            ssht_plot_mollweide(reshape(f_cur{j-J_min+1}(en,:), band_limit , 2*band_limit -1), band_limit , 'Mode', 1);
            %
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Curvelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)], 'FontSize', 8)
        end
	end
end
colormap(jet)
fname = [pltroot,'/s2let_demo8_', configstr, '_curvelet_EarthTopo_MultiRes.png']
print('-r200', '-dpng', fname)


% ---------- 
% Plot reconstructed map (to compare with the initial map): 
% ---------- 
% -----------------
% Signal synthesis: 
% -----------------
f_rec = s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                  'B', B, 'L', L, ...
                                                  'J_min', J_min, ...
                                                  'Spin', Spin, ...
                                                  'Reality', reality, ...
                                                  'Upsample', false, ...
                                                  'SpinLowered', false, ...
                                                  'SpinLoweredFrom', 0);

figure('Position',[100 100 900 200]) 
subplot(2, 2, 1);
% --- plot initial data --- % 
ssht_plot_mollweide(f_gen,L, 'Mode', 1);
title('initial signal')
hold on
subplot(2, 2, 2);
% --- plot reconstructed data --- % 
ssht_plot_mollweide(f_rec,L, 'Mode', 1);
title('reconstructed signal') 
% Check error:
check_error = max(abs(f_gen(:)-f_rec(:)))
                           
fname = [pltroot,'/s2let_demo8_', configstr, '_curvelet_EarthTopo_MultiRes_Int_Rec_signal.png']
print('-r200', '-dpng', fname)

