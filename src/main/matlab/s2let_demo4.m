% s2let_demo4
% Analyse Earth tomography data as a real MW map.
% Compute the wavelet maps and plot them.
% Plot 1 : multiresolution wavelet scales
% Plot 2 : full resolution wavelet scales
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details


load('EGM2008_Topography_flms_L0128');
L = 32;   % To see the multi-resoluition effect i.e. Upsample = false, L >=64) 
flm = flm(1:L^2,1);
Spin = 0; 
f = ssht_inverse(flm, L, 'Reality', true);


%{
inputfile = '../../../data/earth_tomo_mw_128.fits';
[f, L] = s2let_mw_read_real_map(inputfile);
%}

B = 2;
N = L; %to compare with curvelets
J_min = 3; %1
J = s2let_jmax(L, B);

zoomfactor = 1.;  %1.6;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 16; %4 ;  % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 4; % 4 % ns;

maxfigs = nx*ny;
pltroot = '../../../figs'
configstr = ['N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]

[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'B', B, 'J_min', J_min, 'N', N, 'Reality',true, 'Upsample', true, 'Spin', 0);

% For debugging: 
%{
type = 'colour';
zoomfactor = 1.4;
plot_caxis_scale = 1; %2; 
%
figure('Position',[100 100 300 300])
h=subplot(1, 1, 1);
% Plot initial signals on the sphere 
ssht_plot_sphere(f_wav, L, 'Type', type,...
                 'ColourBar', false, 'Lighting', true);
title(h,'Earth map')
locate = get(h,'title');
pos = get(locate,'position');
pos(1,2) = pos(1,2)+0.7;
pos(1,1) = pos(1,1)-0.7;
set(locate,'pos',pos);
zoom(1.2)
v = caxis;
temp = max(abs(v));
caxis([-temp temp]*plot_caxis_scale)
% Plot function on SO(3) by plotting a sphere for each value of gamma
figure('Position',[100 100 1200 1800])
for j = J_min:J 
  % band_limit = min([ s2let_bandlimit(j,J_min,B,L) L ]);  
  %% Compute sampling grids where ngamma = 2*band_limit -1 
  % [alphas, betas, gammas, n, nalpha, nbeta, ngamma] = so3_sampling(band_limit, band_limit, 'Grid', true);
  for i = 1:N    %ngamma,
      h = subplot(ceil(N/6),6,i);
      ssht_plot_sphere(real(f_wav{j-J_min+1,i}), L, 'Type', type,...
                       'ColourBar', false, 'Lighting', true);
      title(h, sprintf('g = %d; Re(f)', i-1))
      locate = get(h,'title');
      pos = get(locate,'Position');
      pos(1,2) = pos(1,2)+0.7;
      pos(1,1) = pos(1,1)-0.7;
      set(locate,'pos',pos);
      v = caxis;
      temp = max(abs(v));
      caxis([-temp temp]*plot_caxis_scale)
      zoom(zoomfactor) 
  end 
end
%}

% FULL-RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(ny, nx, 1);
ssht_plot_mollweide(f, L, 'Mode', 1);
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
fname = [pltroot,'/s2let_demo4_', configstr, '_earth_multires.png']
print('-r200', '-dpng', fname)

% For debugging: 
stop 

[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'B', B, 'J_min', J_min, 'N', N, 'Upsample', false, 'Spin', 0);

% MULTI-RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(ny, nx, 1);
ssht_plot_mollweide(f, L, 'Mode', 1);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(ny, nx, 2);
bl = min([ s2let_bandlimit(J_min-1,J_min,B,L) L ]);
ssht_plot_mollweide(f_scal, bl, 'Mode', 1);
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
            bl =  min([ s2let_bandlimit(j,J_min,B,L) L ]);
            ssht_plot_mollweide(f_wav{j-J_min+1,en}, bl, 'Mode', 1);
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Wavelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)])
        end
	end
end


fname = [pltroot,'/s2let_demo4_', configstr, '_earth_fullres.png']
print('-r200', '-dpng', fname)



