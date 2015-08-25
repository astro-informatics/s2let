% s2let_demo8_Curvelet_Mollweide_randSignals_CompareToOthers
%
% Plot (axiymmetric and directional) wavelet/curvelet coefficients 
% on multiple Mollweide projections.
% The function generates one plot of the scaling function 
% contribution and a grid of plots for each orientation of 
% each scale of the wavelet/curvelet contributions. 
%
% -----------------------------------------------------------
% f_wav or f_cur is cell array with all the wavelet/curvelet coefficients.
% Its first index is the wavelet scale j, the second
% index is the orientation g, and each element is a
% function on the sphere in MW sampling. 
% scal is the corresponding scaling function contribution
% (i.e. just a single function on the sphere).
%
% B is the wavelet parameter.
% L is the angular band-limit.
% N is the orientational band-limit.
% J_min is the first wavelet scale in f_wav/f_cur.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%  'B'               = { Dilation factor; B > 1 (default = 2) }
%  'L'               = { Harmonic band-limit; L > 0 (default = Lguessed) }
%  'J_min'           = { the minimal wavelet scale,(default = 0)}
%  'Spin'            = { Spin number; Spin >= 0 (default = 0) }
%  'Reality'         = { false   [do not assume corresponding signal f real (default)],
%                        true    [assume f real (improves performance)] }
%  'Upsample'        = { false   [multiresolution algorithm (default)],
%                        true    [full resolution wavelets] },
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
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2015  Boris Leistedt, Martin Büttner,
%                     Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

% ---------------
% Define curvelet parameters: 
% ---------------
Spin = 0;             % Spin value of wavelet (for comparison with Axisymmetric wavelets, set to 0)
B = 2;                % B=2 for dyadic sampling
L = 32;               % Angular band-limit (set L>=64 to see multi-resolution effects)
N= L;                 % Azimuthal band-limit (for comparison with curvelets, the azimuthal band limit N = L)
J_min = 3;            % Minimum scale probed by wavelets 
J =s2let_jmax(L, B);  % Maximum scale probed by wavelets =ceil(log L/ log B);

% ---------------
% Generate random complex signals :
% ---------------

disp('Generates random band-limited function')
flm_gen = zeros(L^2,1);
flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');

% -----------------
% Curvelet analysis: (pixel to curvelet space)
% -----------------
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_gen,  ...
                                                          'B', B, 'L', L,  'J_min', J_min, ...
                                                          'Spin', Spin,  ...
                                                          'Reality', false,...
                                                          'Upsample', true);

%
% FULL RESOLUTION PLOT (Upsample is true)
%
zoomfactor = 1.;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 10; % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 4;  % ns;

maxfigs = nx*ny;
pltroot = '../../../figs' ; 
configstr = ['N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]; 

figure('Position',[100 100 1300 1000])
subplot(ny, nx, 1);
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(ny, nx, 2);
%
ssht_plot_mollweide(f_scal, L, 'Mode', 1);
%
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
title('Scaling fct')
ind = 2; 
for j = J_min:J
	for en = 1:N-1
		ind = ind + 1 ; 
        if ind <= maxfigs
            subplot(ny, nx, ind);
            ssht_plot_mollweide(reshape(f_cur{j-J_min+1}(en,:), L, 2*L-1), L, 'Mode', 1);
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Wavelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)])
        end
	end
end
colormap(jet)
fname = [pltroot,'/s2let_demo8_', configstr, '_curvelet_Mollweide_randSignals.png']
print('-r200', '-dpng', fname)



% ============
% Directional wavelet analysis (pixel to wavelet space)
% ============
[f_wav, f_scal] = s2let_transform_analysis_mw(f_gen, 'B', B, 'J_min', J_min, ...
                                              'N', N, 'Reality',false,...
                                              'Upsample', true, 'Spin', 0);
%
% FULL RESOLUTION PLOT (Upsample is true)
%
figure('Position',[100 100 1300 1000])
subplot(ny, nx, 1);
%
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
%
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(ny, nx, 2);
%
ssht_plot_mollweide(f_scal, L, 'Mode', 1);
%
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
title('Scaling fct')
ind = 2;
for j = J_min:J
	for en = 1:N-1
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
fname = [pltroot,'/s2let_demo8_', configstr, '_directional_wavelet_Mollweide_randSignals.png']
print('-r200', '-dpng', fname)


zoomfactor = 1.;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = ns;
% ============
% Axisymmetric wavelet analysis  (pixel to wavelet space)
% ============
[f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f_gen, 'B', B, 'J_min', J_min, ...
                                                     'Reality', false, 'Upsample', true);
%
% FULL RESOLUTION PLOT (Upsample is true)
%
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
%
ssht_plot_mollweide(f_gen, L);
%
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(nx, ny, 2);
%
ssht_plot_mollweide(f_scal, L);
%
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
title('Scaling fct')
for j = J_min:J
   subplot(nx, ny, j-J_min+3);
   %
   ssht_plot_mollweide(f_wav{j-J_min+1}, L);
   %
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
   %title(['Wavelet scale : ',int2str(j)-J_min+1])
end
colormap(jet)
fname = [pltroot,'/s2let_demo8_', configstr, '_axisymmetric_wavelet_Mollweide_randSignals.png']
print('-r200', '-dpng', fname)

