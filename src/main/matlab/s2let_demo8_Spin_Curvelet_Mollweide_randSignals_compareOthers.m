% s2let_demo8_Curvelet_Mollweide_randSignals
% -----------------------------------------------------------
% Plot curvelet coefficients on multiple Mollweide projections.
% The function generates one plot of the scaling function
% contribution and a grid of plots for each orientation of
% each scale of the wavelet contributions. 
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
%
%  'Upsample'      = { false   [multiresolution algorithm (default)],
%                      true  [full resolution wavelets] },
%  'Function'        = { 'real' [plot the real part of the input functions (default)],
%                        'imag' [plot the imaginary part of the input functions],
%                        'abs'  [plot the absolute value of the input functions] }
% -----------------------------------------------------------
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
%
% Modified S2LET package to perform Curvelets on the Sphere.
% -----------------------------------------------------------

% ---------------
% Define curvelet parameters: 
% ---------------
Spin = 0;
B = 2;   % for dyadic sampling 
L = 32;   % To see the multi-resoluition effect i.e. Upsample = false, L >=64) 
N= L;     % Since m=l, the azimuthal band limit N = overall band limit L
J_min = 3; % minimum and maximum scale probed by wavelets 
J =s2let_jmax(L, B);  %=ceil(log L/ log B);  

% ---------------
% Generate random complex signals :
% ---------------

disp('Generates random band-limited function')
flm_gen = zeros(L^2,1);
flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');
%{
load('EGM2008_Topography_flms_L0128');
flm_gen = flm(1:L^2,1);
f_gen = ssht_inverse(flm_gen, L,'Reality', true);
%}

%{
% -----------------
% Signal analysis: (harmonic to curvelet space) 
% -----------------
% Call matlab function analysis_lm2cur
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_lm2cur(flm_gen,  ...
                                                          'B', B, 'L', L,  'J_min', J_min, ...
                                                          'Spin', Spin,  ...
                                                          'Sampling', 'MW',...
                                                          'Upsample', true);
                                                     
f_cur_new = cell(J+1-J_min, N-1);  
for j = J_min:J 
   for en = 1: N-1 
   f_cur_new{j-J_min+1, en} = reshape(f_cur{j-J_min+1}(en,:), L, 2*L-1);
   end
end   

%
% FULL RESOLUTION PLOT
%
zoomfactor = 1.;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 10; % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 4; % ns;

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
            ssht_plot_mollweide(f_cur_new{j-J_min+1,en}, L, 'Mode', 1);
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Wavelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)])
        end
	end
end
colormap(jet)
%fname = [pltroot,'/s2let_demo9_', configstr, '_testSIGANL.png']
%print('-r200', '-dpng', fname)
%}
% -----------------
% Signal analysis: (pixel to curvelet space) 
% -----------------
% Call matlab function analysis_lm2cur
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_px2cur(f_gen,  ...
                                                          'B', B, 'L', L,  'J_min', J_min, ...
                                                          'Spin', Spin,  ...
                                                          'Reality', false,...
                                                          'Sampling', 'MW',...
                                                          'Upsample', true);

%
% FULL RESOLUTION PLOT
%
zoomfactor = 1.;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 10; % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 4; % ns;

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
fname = [pltroot,'/s2let_demo8_', configstr, '_curvelets.png']
print('-r200', '-dpng', fname)



% ============
% Directional wavelets
% ============

[f_wav, f_scal] = s2let_transform_analysis_mw(f_gen, 'B', B, 'J_min', J_min, ...
                                              'N', N, 'Reality',false,...
                                              'Upsample', true, 'Spin', 0);

% FULL-RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(ny, nx, 1);
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(ny, nx, 2);
ssht_plot_mollweide(f_scal, L, 'Mode', 1);
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
fname = [pltroot,'/s2let_demo8_', configstr, '_directional_wavelets.png']
print('-r200', '-dpng', fname)



% ============
% Axisymmetric wavelets
% ============
zoomfactor = 1.;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = ns;
% Perform decomposition
[f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f_gen, 'B', B, 'J_min', J_min, ...
                                                     'Reality', false, 'Upsample', true);
% FULL RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
ssht_plot_mollweide(f_gen, L);
%title('Initial data')
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
fname = [pltroot,'/s2let_demo8_', configstr, '_axisymmetric_wavelets.png']
print('-r200', '-dpng', fname)

