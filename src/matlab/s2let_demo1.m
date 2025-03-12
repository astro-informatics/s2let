% s2let_demo1
% Analyse Earth tomography data as a real MW map.
% Compute the wavelet maps and plot them.
% Plot 1 : multiresolution wavelet scales
% Plot 2 : full resolution wavelet scales
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

load('EGM2008_Topography_flms_L0128');
f = ssht_inverse(flm, L, 'Reality', true);

%inputfile = 'data/earth_tomo_mw_128.fits';
%[f, L] = s2let_mw_read_real_map(inputfile);

B = 3;
J_min = 2;
J = s2let_jmax(L, B);

zoomfactor = 1.2;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = ns;

% Perform decomposition
[f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f, 'B', B, 'J_min', J_min, 'Reality', true, 'Upsample', true);
% FULL RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
ssht_plot_mollweide(f, L);
%title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(nx, ny, 2);
ssht_plot_mollweide(f_scal, L);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
%title('Scaling fct')
for j = J_min:J
   subplot(nx, ny, j-J_min+3);
   ssht_plot_mollweide(f_wav{j-J_min+1}, L);
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
   %title(['Wavelet scale : ',int2str(j)-J_min+1])
end


% Perform decomposition
[f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f, 'B', B, 'J_min', J_min, 'Reality', true);

% MULTIRESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
ssht_plot_mollweide(f, L);
%title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(nx, ny, 2);
bl = min([ s2let_bandlimit(J_min-1,J_min,B,L) L ]);
ssht_plot_mollweide(f_scal, bl);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
%title('Scaling fct')
for j = J_min:J
   subplot(nx, ny, j-J_min+3);
   bl =  min([ s2let_bandlimit(j,J_min,B,L) L ]);
   ssht_plot_mollweide(f_wav{j-J_min+1}, bl);
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
   %title(['Wavelet scale : ',int2str(j)-J_min+1])
end
