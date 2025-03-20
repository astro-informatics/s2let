% s2let_demo4
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

B = 6;
N = 3;
J_min = 1;
J = s2let_jmax(L, B);

zoomfactor = 1.8;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 4; % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 3; % ns;

maxfigs = nx*ny;
pltroot = '../../../figs'
configstr = ['N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]

[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'B', B, 'J_min', J_min, 'N', N, 'Upsample', true, 'Spin', 0);

% FULL RESOLUTION PLOT
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
fname = [pltroot,'/s2let_demo4_', configstr, '_earth_fullres.png']
print('-r200', '-dpng', fname)




[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'B', B, 'J_min', J_min, 'N', N, 'Upsample', false, 'Spin', 0);

% FULL RESOLUTION PLOT
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


fname = [pltroot,'/s2let_demo4_', configstr, '_earth_multires.png']
print('-r200', '-dpng', fname)

