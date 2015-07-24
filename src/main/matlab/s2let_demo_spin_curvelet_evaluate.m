
% s2let_demo_spin_curvelet_evaluate
% Evaluate timing and error of spin curvelet transform.
%
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear;

N_test = 3 %10

B = 2;
J_min = 0;
spin = 0;
reality = false;
sampling_method = 'MW';
save_plots = true;

Ls = [4 8 16 ]   %[32 64 128 256 512]

err = zeros(N_test, length(Ls));
time_analysis = zeros(N_test, length(Ls));
time_synthesis = zeros(N_test, length(Ls));

el_ind = 0;
for L = Ls
   el_ind = el_ind + 1
   
   for n = 1:N_test
   
      n
      
      %disp('Generate band-limited function in harmonic space')
      flm = zeros(L^2,1);
      flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
      flm = 2.*(flm - (1+sqrt(-1))./2);

      %disp('Construct the corresponding signal on the sphere')
      f = ssht_inverse(flm, L, 'Method', sampling_method);
      
      tstart = tic;
      [f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_px2cur(f, ...
         'B', B, 'L', L, 'J_min', J_min, ...
         'Upsample', true, 'Spin', spin, 'Reality', reality, ...
         'Sampling', sampling_method);
      time_analysis(n, el_ind) = toc(tstart);
      
      tstart = tic;
      f_recov = s2let_transform_spin_curvelet_synthesis_cur2px(f_cur, f_scal, ...
         'L', L, 'B', B, 'J_min', J_min, ...
         'Upsample', true, 'Spin', spin, 'Reality', reality, ...
         'Sampling', sampling_method);


      time_synthesis(n, el_ind) = toc(tstart);

      err(n, el_ind) = max(abs(f(:) - f_recov(:)));

   end

end

err_mean = mean(err);
err_std = std(err);
err_std_log = std(log10(err));

time_analysis_mean = mean(time_analysis);
time_analysis_std = std(time_analysis);


time_synthesis_mean = mean(time_synthesis);
time_synthesis_std = std(time_synthesis);


time_total = time_analysis + time_synthesis;
time_total_mean = mean(time_total);
time_total_std = std(time_total);
time_total_std_log = std(log10(time_total));


%% Define saved figure location 
pltroot = '../../../figs/' ;

%% Define plotting parameters.

istart = 1;
iend = length(Ls);

line_width = 1.8;
line_width_thick = 2.5;
marker_size = 7;
marker_type = 'o';
green_light = [0.2 0.6 0.4];% x339966
green_dark = [0 0.4 0.2];   % x006633
blue_light = [0.2 0.4 0.8]; % x3366CC
blue_dark = [0 0 1];        % x0000FF
red_light = [1 0.4 0.2];    % xFF6633
red_dark = [0.8 0.2 0];     % xCC3300


% Plot error
figure;
   
plot(log2(Ls(istart:iend)), ...
   log10(err_mean(istart:iend)), ...
   'Color', green_dark, ...
   'Marker', marker_type, ...
   'MarkerSize', marker_size, ...
   'MarkerFaceColor', green_light, ...
   'MarkerEdgeColor', green_dark, ...
   'LineStyle', '--', ...
   'LineWidth', line_width);

hold on;

plot(log2(Ls(istart:iend)), ...
   log10(Ls(istart:iend).^2)-15, ...
   'r', ...
   'LineWidth', line_width_thick);
a = axis;
set(gca,'XTick',a(1):1:a(2));
set(gca,'XTickLabel', 2.^[a(1):1:a(2)]);
set(gca,'YTick',floor(a(3)):1:ceil(a(4)));
set(gca,'YTickLabel',{10.^[floor(a(3)):ceil(a(4))]});
xlabel('L');
ylabel('E');

set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 20)

% if save_plots, print('-depsc2', 'plots/spin_curvelet_error.eps'); end

fname = [pltroot,'spin_curvelet_error.png']
print('-r200', '-dpng', fname);

% Plot error with error bars
figure
  
errorbar(log2(Ls(istart:iend)), ...
   log10(err_mean(istart:iend)), ...
   err_std_log(istart:iend), ...
   'Color', green_dark, ...
   'Marker', marker_type, ...
   'MarkerSize', marker_size, ...
   'MarkerFaceColor', green_light, ...
   'MarkerEdgeColor', green_dark, ...
   'LineStyle', '--', ...
   'LineWidth', line_width);

hold on;

plot(log2(Ls(istart:iend)), ...
   log10(Ls(istart:iend).^2)-15, ...
   'r', ...
   'LineWidth', line_width_thick);
a = axis;
set(gca,'XTick',a(1):1:a(2));
set(gca,'XTickLabel', 2.^[a(1):1:a(2)]);
set(gca,'YTick',floor(a(3)):1:ceil(a(4)));
set(gca,'YTickLabel',{10.^[floor(a(3)):ceil(a(4))]});
xlabel('L');
ylabel('E');

set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 20)

% if save_plots, print('-depsc2', 'plots/spin_curvelet_error_errorbars.eps'); end
fname = [pltroot,'spin_curvelet_errorbar.png']
print('-r200', '-dpng', fname);

% Plot timing
figure;
   
plot(log2(Ls(istart:iend)), ...
   log10(time_total_mean(istart:iend)), ...
   'Color', green_dark, ...
   'Marker', marker_type, ...
   'MarkerSize', marker_size, ...
   'MarkerFaceColor', green_light, ...
   'MarkerEdgeColor', green_dark, ...
   'LineStyle', '--', ...
   'LineWidth', line_width);

hold on;

plot(log2(Ls(istart:iend)), ...
   log10(Ls(istart:iend).^3)-5, ...
   'r', ...
   'LineWidth', line_width_thick);
a = axis;
set(gca,'XTick',a(1):1:a(2));
set(gca,'XTickLabel', 2.^[a(1):1:a(2)]);
set(gca,'YTick',floor(a(3)):1:ceil(a(4)));
set(gca,'YTickLabel',{10.^[floor(a(3)):ceil(a(4))]});
xlabel('L');
ylabel('T');

set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 20)

% if save_plots, print('-depsc2', 'plots/spin_curvelet_timing.eps'); end

fname = [pltroot,'spin_curvelet_timing.png']
print('-r200', '-dpng', fname);

% Plot timing with error bars
figure;
   
errorbar(log2(Ls(istart:iend)), ...
   log10(time_total_mean(istart:iend)), ...
   time_total_std_log(istart:iend), ...
   'Color', green_dark, ...
   'Marker', marker_type, ...
   'MarkerSize', marker_size, ...
   'MarkerFaceColor', green_light, ...
   'MarkerEdgeColor', green_dark, ...
   'LineStyle', '--', ...
   'LineWidth', line_width);

hold on;

plot(log2(Ls(istart:iend)), ...
   log10(Ls(istart:iend).^3)-5, ...
   'r', ...
   'LineWidth', line_width_thick);
a = axis;
set(gca,'XTick',a(1):1:a(2));
set(gca,'XTickLabel', 2.^[a(1):1:a(2)]);
set(gca,'YTick',floor(a(3)):1:ceil(a(4)));
set(gca,'YTickLabel',{10.^[floor(a(3)):ceil(a(4))]});
xlabel('L');
ylabel('T');

set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 20)

% if save_plots, print('-depsc2', 'plots/spin_curvelet_timing_errorbars.eps'); end

fname = [pltroot,'spin_curvelet_timing_errorbar.png']
print('-r200', '-dpng', fname);

