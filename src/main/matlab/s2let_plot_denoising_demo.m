%% Preparation

% How many sigmas of the spin signal data should be contained in the
% color bar.
sigmas = 2;

position = [ 500 500 1000 800 ];
set(0, 'DefaultFigurePosition', position);
plot_root = '../../../figs/';
plot_size = [0 0 23 13];

% Set working directory to script's location
% (assumed to be src/main/matlab)

cd(fileparts(mfilename('fullpath')))

%% Plot spin map results

[f_Q, L] = s2let_mw_read_real_map('../../../data/spin_signal_real_input.fits');
[f_Q_noise, ~] = s2let_mw_read_real_map('../../../data/spin_signal_real_noise.fits');
[f_Q_input_noise, ~] = s2let_mw_read_real_map('../../../data/spin_signal_real_input_noise.fits');
[f_Q_denoised, ~] = s2let_mw_read_real_map('../../../data/spin_signal_real_denoised.fits');

f_error = log(abs(f_Q_denoised./f_Q - 1));

mean_f = mean(f_Q_input_noise(:));
stddev_f = std(f_Q_input_noise(:));
range = [mean_f-sigmas*stddev_f mean_f+sigmas*stddev_f];

figure
ssht_plot_mollweide(f_Q, L, 'Mode', 1)
colorbar;
caxis(range);
title('Spin signal (Q) - Input map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_Q_noise, L, 'Mode', 1)
colorbar;
caxis(range);
title('Spin signal (Q) - Noise map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_noise.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_Q_input_noise, L, 'Mode', 1)
colorbar;
caxis(range);
title('Spin signal (Q) - Input map with added noise');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_input_noise.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_Q_denoised, L, 'Mode', 1)
colorbar;
caxis(range);
title('Spin signal (Q) - Denoised map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
colorbar;
title('Spin signal (Q) - Logarithmic relative error');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_error.png'];
print('-r200', '-dpng', fname)

[f_U, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_input.fits');
[f_U_noise, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_noise.fits');
[f_U_input_noise, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_input_noise.fits');
[f_U_denoised, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_denoised.fits');

f_error = log(abs(f_U_denoised./f_U - 1));

mean_f = mean(f_U_input_noise(:));
stddev_f = std(f_U_input_noise(:));
range = [mean_f-sigmas*stddev_f mean_f+sigmas*stddev_f];

figure
ssht_plot_mollweide(f_U, L, 'Mode', 1)
colorbar;
caxis(range);
title('Spin signal (U) - Input map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_U_noise, L, 'Mode', 1)
colorbar;
caxis(range);
title('Spin signal (U) - Noise map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_noise.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_U_input_noise, L, 'Mode', 1)
colorbar;
caxis(range);
title('Spin signal (U) - Input map with added noise');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_input_noise.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_U_denoised, L, 'Mode', 1)
colorbar;
caxis(range);
title('Spin signal (U) - Denoised map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
colorbar;
title('Spin signal (U) - Logarithmic relative error');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_error.png'];
print('-r200', '-dpng', fname)

%% Plot Q + iU map

subL = 64;

f_QU = f_Q + 1i*f_U;
f_QU_noise = f_Q_noise + 1i*f_U_noise;
f_QU_input_noise = f_Q_input_noise + 1i*f_U_input_noise;
f_QU_denoised = f_Q_denoised + 1i*f_U_denoised;

f_error = log(abs(f_QU_denoised./f_QU - 1));

mean_f = abs(mean(f_QU_input_noise(:)));
stddev_f = abs(std(f_QU_input_noise(:)));
range = [0 mean_f+sigmas*stddev_f];

figure
ssht_plot_mollweide(f_QU, L, 'Mode', 3, 'SubL', subL)
colorbar;
caxis(range);
title('Spin signal (Q + iU) - Input map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_QU_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_QU_noise, L, 'Mode', 3, 'SubL', subL)
colorbar;
caxis(range);
title('Spin signal (Q + iU) - Noise map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_QU_noise.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_QU_input_noise, L, 'Mode', 3, 'SubL', subL)
colorbar;
caxis(range);
title('Spin signal (Q + iU) - Input map with added noise');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_QU_input_noise.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_QU_denoised, L, 'Mode', 3, 'SubL', subL)
colorbar;
caxis(range);
title('Spin signal (Q + iU) - Denoised map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_QU_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
colorbar;
title('Spin signal (Q + iU) - Logarithmic relative error');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_QU_error.png'];
print('-r200', '-dpng', fname)


%% Plot real map results

[f, L] = s2let_mw_read_real_map('../../../data/real_signal_input.fits');
[f_noise, ~] = s2let_mw_read_real_map('../../../data/real_signal_noise.fits');
[f_input_noise, ~] = s2let_mw_read_real_map('../../../data/real_signal_input_noise.fits');
[f_denoised, ~] = s2let_mw_read_real_map('../../../data/real_signal_denoised.fits');

% Precompute Wigner small-d functions
sqrt_tbl = sqrt([0:2*(L-1)+1])';
signs = ones(L+1,1);
signs(2:2:end) = -1;

d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, 0, ...
  'SqrtTable', sqrt_tbl, 'SignTable', signs);
for el = 1:L-1
  d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, 0, ...
     'SqrtTable', sqrt_tbl, 'SignTable', signs);
end

% Rotate maps by 180°
f = ssht_inverse(ssht_rotate_flm(ssht_forward(f, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_noise = ssht_inverse(ssht_rotate_flm(ssht_forward(f_noise, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_input_noise = ssht_inverse(ssht_rotate_flm(ssht_forward(f_input_noise, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_denoised = ssht_inverse(ssht_rotate_flm(ssht_forward(f_denoised, L, 'Reality', true), d, pi, 0), L, 'Reality', true);

f_error = log(abs(f_denoised./f - 1));

figure
ssht_plot_mollweide(f, L, 'Mode', 1)
colorbar;
title('Real signal - Input map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_noise, L, 'Mode', 1)
colorbar;
title('Real signal - Noise map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_noise.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_input_noise, L, 'Mode', 1)
colorbar;
title('Real signal - Input map with added noise');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_input_noise.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_denoised, L, 'Mode', 1)
colorbar;
title('Real signal - Denoised map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
colorbar;
title('Real signal - Logarithmic relative error');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_error.png'];
print('-r200', '-dpng', fname)


%% Clean up

close all

clear position plot_root plot_size
clear f f_noise f_input_noise f_denoised f_error
clear f_Q f_Q_noise f_Q_input_noise f_Q_denoised
clear f_U f_U_noise f_U_input_noise f_U_denoised
clear f_QU f_QU_noise f_QU_input_noise f_QU_denoised
clear L
clear fname
clear d el signs sqrt_tbl
clear sigmas mean_f stddev_f range
