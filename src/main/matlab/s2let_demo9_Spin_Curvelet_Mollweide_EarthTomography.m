% s2let_demo9_Spin_Curvelet_Mollweide_EarthTomography
% -----------------------------------------------------------
% Plot curvelet coefficients on multiple Mollweide projections.
% The function generates one plot of the scaling function
% contribution and a grid of plots for each orientation of
% each scale of the wavelet contributions. 
% -----------------------------------------------------------
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% 
% Modified S2LET package to perform Curvelets on the Sphere.
% -----------------------------------------------------------
% wav is cell array with all the wavelet coefficients.
% its first index is the wavelet scale j, the second
% index is the orientation g, and each element is a
% function on the sphere in MW sampling.
% scal is the corresponding scaling function contribution
% (i.e. just a single function on the sphere).
% B is the wavelet parameter.
% L is the angular band-limit.
% N is the orientational band-limit.
% J_min is the first wavelet scale in wav.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%  'Spin'            = { non-negative integers (default=0) }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'Upsample'        = { false        [multiresolution algorithm (default)],
%                        true       [full resolution curvelets] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               curvelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the curvelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the curvelets
%                       should be lowered from (default = 0)]
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%
%
% FOR ssht_plot_mollweide(f, L, <options>)
% where f is the sampled function and L is the harmonic band-limit.
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Method'          = { 'MW'         [McEwen & Wiaux sampling (default)],
%                        'MWSS'       [McEwen & Wiaux symmetric sampling],
%                        'DH'         [Driscoll & Healy sampling],
%                        'GL'         [Gauss-Legendre sampling] }
%  'ColourBar'       = { false        [do not add colour bar (default)],
%                        true         [add colour bar] }
%  'Mode'            = { 0            Plot amplitude, or modulus is f complex (default),
%                        1            Plot real part,
%                        2            Plot imaginaty part,
%                        3            Plot modulus and arrows for real/img angle }
%  'Spin'            = { non-negative integers (default=0) }
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2015  
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

clear all;
% close all ;

load('EGM2008_Topography_flms_L0128');
L = 32;   % To see the multi-resoluition effect i.e. Upsample = false, L >=64) 
flm_gen = flm(1:L^2,1);
f_gen = ssht_inverse(flm_gen, L,'Reality', true);

% ---------------
% Define curvelet parameters: 
% ---------------
B = 2;   
Spin = 0; 
N= L;     % Since m=l, the azimuthal band limit N = overall band limit L
J_min = 3; % minimum and maximum scale probed by wavelets 
J =s2let_jmax(L, B);  %=ceil(log L/ log B);  

% ---------------
% Define Plotting  parameters: 
% ---------------
zoomfactor =1.2;  %1.6;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 15; %10  %16 % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 5;  %5 % ns;

maxfigs = nx*ny;
pltroot = '../../../figs' ; 
configstr = ['N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]; 

% ================================================
% FULL RESOLUTION PLOT (Upsample: true)
% ================================================
% Signal analysis: 
% -----------------
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_px2cur(f_gen,  ...
                                                  'B', B, 'L', L, ...
                                                  'J_min', J_min, ...
                                                  'Spin', Spin, ...
                                                  'Reality', true, ...
                                                  'Upsample', true, ...
                                                  'SpinLowered', false, ...
                                                  'SpinLoweredFrom', 0,...
                                                  'Sampling', 'MW');
                                              
f_cur_new = cell(J+1-J_min, N-1);  
for j = J_min:J 
   for en = 1: N-1   %2*N-1 
   f_cur_new{j-J_min+1, en} = reshape(f_cur{j-J_min+1}(en,:), L, 2*L-1);
   end
end   
% figure
% ssht_plot_mollweide(f_cur_new{1, 1}, L, 'Mode', 1);


figure('Position',[20 20 1700 1400]) %100 100 1300 1000
subplot(ny, nx, 1);
% --- plot initial data --- % 
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
%
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
% 
subplot(ny, nx, 2);
% --- plot scaling function contributions --- % 
ssht_plot_mollweide(f_scal, L, 'Mode', 1);
%
title('Scaling fct')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
% --- plot curvelet kernel contributions --- % 
ind = 2;
for j = J_min:J
	for en = 1: N-1 %2*N-1 %N
		ind = ind + 1;
        if ind <= maxfigs
            subplot(ny, nx, ind);
            %
            ssht_plot_mollweide(f_cur_new{j-J_min+1, en}, L, 'Mode', 1);
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
fname = [pltroot,'/s2let_demo9_', configstr, '_spin_curvelet_EarthTomo_fullres.png']
print('-r200', '-dpng', fname)


% ---------- 
% Compare reconstructed signal with the initial signals: 
% ---------- 
f_rec = s2let_transform_spin_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                      'B', B, 'L', L, ...
                                                      'J_min', J_min, ...
                                                      'Spin', Spin, ...
                                                      'Reality', true, ...
                                                      'Upsample', true, ...
                                                      'SpinLowered', false, ...
                                                      'SpinLoweredFrom', 0,...
                                                 	  'Sampling', 'MW');

figure('Position',[100 100 900 200]) 
subplot(2, 2, 1);
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
title('initial signal')
hold on
subplot(2, 2, 2);
ssht_plot_mollweide(f_rec,L, 'Mode', 1);
title('reconstructed signal')
fname = [pltroot,'/s2let_demo9_', configstr, '_spin_curvelet_EarthTomo_Fullres_Int_Rec_signal.png']
print('-r200', '-dpng', fname)
% Check error:
check_error = max(abs(f_gen(:)-f_rec(:)))
                           

% ================================================
% MULTI-RESOLUTION PLOT (Upsample: false)
% ================================================
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_px2cur(f_gen,  ...
                                                  'B', B, 'L', L, ...
                                                  'J_min', J_min, ...
                                                  'Spin', Spin, ...
                                                  'Reality', true, ...
                                                  'Upsample', false, ...
                                                  'SpinLowered', false, ...
                                                  'SpinLoweredFrom', 0,...
                                                  'Sampling', 'MW');  

                                              
figure('Position',[20 20 1700 1400]) %100 100 1300 1000
subplot(ny, nx, 1);
% --- plot initial data --- % 
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
%
title('Initial data')
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
title('Scaling fct')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
% --- plot curvelet kernel contributions --- % 
ind = 2;
for j = J_min:J
	for en = 1: N-1   %2*N-1 %N
		ind = ind + 1;
        if ind <= maxfigs
            subplot(ny, nx, ind);
            %
            bl =  min([ s2let_bandlimit(j,J_min,B,L) L ]);
            ssht_plot_mollweide(reshape(f_cur{j-J_min+1}(en,:), bl, 2*bl-1), bl, 'Mode', 1);
            %
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Curvelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)], 'FontSize', 10)
        end
	end
end
colormap(jet)
fname = [pltroot,'/s2let_demo9_', configstr, '_spin_curvelet_EarthTomo_multires.png']
print('-r200', '-dpng', fname)


% ---------- 
% Compare reconstructed signal with the initial signals: 
% ---------- 
f_rec = s2let_transform_spin_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                       'B', B, 'L', L, ...
                                                       'J_min', J_min, ...
                                                       'Spin', Spin, ...
                                                       'Reality', true, ...
                                                       'Upsample', false, ...
                                                       'SpinLowered', false, ...
                                                       'SpinLoweredFrom', 0,...
                                                       'Sampling', 'MW');

figure('Position',[100 100 900 200]) 
subplot(2, 2, 1);
ssht_plot_mollweide(f_gen,L, 'Mode', 1);
title('initial signal')
hold on
subplot(2, 2, 2);
ssht_plot_mollweide(f_rec,L, 'Mode', 1);
title('reconstructed signal')
% Check error:
check_error = max(abs(f_gen(:)-f_rec(:)))
                           
fname = [pltroot,'/s2let_demo9_', configstr, '_spin_curvelet_EarthTomo_multires_Int_Rec_signal.png']
print('-r200', '-dpng', fname)

