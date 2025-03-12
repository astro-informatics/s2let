
% s2let_demo_ridgelet_evaluate
% Evaluate timing and error of ridgelet transform.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear;

N_test = 3

B = 2;
J_min = 0;
spin = 2;
reality = false;
sampling_method = 'MWSS';
save_plots = false;

Ls = [32 64] % 128 256 512]

err = zeros(N_test, length(Ls));
time_analysis = zeros(N_test, length(Ls));
time_synthesis = zeros(N_test, length(Ls));

el_ind = 0;
for L = Ls
   el_ind = el_ind + 1
   
   for n = 1:N_test
   
      n
      
      % Generate band-limited function in harmonic space
      flm = zeros(L^2,1);
      flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
      flm = 2.*(flm - (1+sqrt(-1))./2);
      
      % For spin signals, set coefficients for el < spin to zero.
      if spin > 0
         for el = 0:abs(spin)-1
            for m = -el:el
               ind = ssht_elm2ind(el, m);
               flm(ind) = 0;
            end
         end
      end
      
      % Impose condition to ensure invertible
      % (imposes antipodal symmetry for scalar signals).
      for el = max([0 abs(spin)]):L-1         
         if mod(el + spin, 2) == 1
            for m = -el:el
               ind = ssht_elm2ind(el, m);
               flm(ind) = 0;
            end
         end
      end
            
      % Impose antipodal symmetry for spin signals.
%       for el = 0:L-1
%          ind = el*el + el + 1;         
%          if mod(el + spin, 2) == 1
%             flm(ind,1) = 0;
%          end
%          for m = 1:el
%             ind_pm = el*el + el + m + 1;
%             ind_nm = el*el + el - m + 1;
%             flm(ind_nm,1) = (-1).^el .* (-1)^m .* (-1)^spin * conj(flm(ind_pm,1));
%          end
%       end                  
      
      % Construct the corresponding signal on the sphere
      f = ssht_inverse(flm, L, 'Method', sampling_method, 'Spin', spin);                  
      
      tstart = tic;
      [f_ridgelet_wav, f_ridgelet_scal] = s2let_ridgelet_analysis(f, ...
         'B', B, 'J_min', J_min, ...
         'Upsample', true, 'Spin', spin, 'Reality', reality, ...
         'Sampling', sampling_method);
      time_analysis(n, el_ind) = toc(tstart);
      
      tstart = tic;
      f_recov = s2let_ridgelet_synthesis(f_ridgelet_wav, f_ridgelet_scal, ...
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

if save_plots, print('-depsc2', 'plots/ridgelet_error.eps'); end

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

if save_plots, print('-depsc2', 'plots/ridgelet_error_errorbars.eps'); end


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

if save_plots, print('-depsc2', 'plots/ridgelet_timing.eps'); end


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

if save_plots, print('-depsc2', 'plots/ridgelet_timing_errorbars.eps'); end


