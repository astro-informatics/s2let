
mapTfile = '../../../../ebsep/data/wmap_T_rot.fits';
mapQfile = '../../../../ebsep/data/wmap_Q_rot.fits';
mapUfile = '../../../../ebsep/data/wmap_U_rot.fits';
maskTfile = '../../../../ebsep/data/wmap_T_mask_rot.fits';
maskQUfile = '../../../../ebsep/data/wmap_QU_mask_rot.fits';

[mapT_hpx, ~] = s2let_hpx_read_real_map(mapTfile);
[maskT_hpx, ~] = s2let_hpx_read_real_map(maskTfile);
[mapQ_hpx, ~] = s2let_hpx_read_real_map(mapQfile);
[mapU_hpx, ~] = s2let_hpx_read_real_map(mapUfile);
[maskQU_hpx, ~] = s2let_hpx_read_real_map(maskQUfile);

Spin = 2 
L = 192
B = 3
N = 3
J_min = 2
downsampling = false;
mode = 3;
PhiRot = pi;
threshold = -1;

% Analyse Q/U maps with Healpix routines to produce full-sky E/B alms
[E_lm, B_lm] = s2let_hpx_map2alm_spin(mapQ_hpx, mapU_hpx, 'L', L);

% Compute full sky Q+iU and Q-iU maps from E/B alms
[QpU_lm] = ebsep_compute_QpU_fullsky(E_lm, B_lm);
[QmU_lm] = ebsep_compute_QmU_fullsky(E_lm, B_lm);
QpU = ssht_inverse(QpU_lm, L, 'spin', 2);
QmU = ssht_inverse(QmU_lm, L, 'spin', -2);
maskQU = s2let_hpx2mw(maskQU_hpx, 'L', L);
maskQU(maskQU > threshold) = 1;
maskQU(maskQU <= threshold) = 0;
maskQU_lm = ssht_forward(maskQU, L);

[maskQU, mask_wav, mask_scal] = ebsep_sample_mask(maskQU_lm, B, L, J_min, 'threshold', threshold, 'Downsample', downsampling);
%
f = QpU .* maskQU;
%f = QpU;
%

J = s2let_jmax(L, B);
zoomfactor = 1.7;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 3; % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 3; % ns;

maxfigs = nx*ny;
pltroot = '../../../figs'
configstr = ['Spin',int2str(Spin),'_N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]

[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'B', B, 'J_min', J_min, 'N', N, 'Upsample', ~downsampling);

bl = L;
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
ssht_plot_mollweide(f, L, 'Mode', mode);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(nx, ny, 2);
if downsampling
    bl = min([ s2let_bandlimit(J_min-1,J_min,B,L), L ]);
end
ssht_plot_mollweide(f_scal .* mask_scal, bl, 'Mode', mode);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
title('Scaling fct')
ind = 3
for j = J_min:J
	for en = 1:N
		ind = ind + 1
        if ind <= maxfigs
            subplot(nx, ny, ind);
            if downsampling
                bl = min([s2let_bandlimit(j, J_min, B, L), L]);
            end
            ssht_plot_mollweide(f_wav{j-J_min+1,en}.*mask_wav{j-J_min+1,1}, bl, 'Mode', mode);
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Wavelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)])
         end
	end
end

fname = [pltroot,'/s2let_demo6_', configstr, '_wmap_unmasked.png']
print('-r200', '-dpng', fname)


