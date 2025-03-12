
% s2let_demo_ridgelet_plot
% Compute and plot the ridgelet wavelet and scaling functions.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

% close all;

L = 256
B = 2;
J_min = 3;
J = s2let_jmax(L, B);
sampling_method = 'MWSS';
% sphere_plot_type = 'parametric';
sphere_plot_type = 'colour';
sphere_plot_lighting = true;
sphere_plot_scale = 1;
save_plots = false;
spin = 2
reality = false

% plot_component = @(x) (x);
plot_component = @(x) real(x);
% plot_component = @(x) imag(x);
% plot_component = @(x) abs(x);

[ridgelet_wav, ridgelet_scal] = s2let_ridgelet_compute_wav(L, ...
   'B', B, 'J_min', J_min, ...
   'Spin', spin, 'Reality', reality, 'Sampling', sampling_method);

figure
ssht_plot_sphere(plot_component(ridgelet_scal), L, ...
   'Lighting', sphere_plot_lighting, ...
   'Type', sphere_plot_type, 'ParametricScale', sphere_plot_scale, ...
   'Method', sampling_method, ...
   'ColourBar', true);
c = caxis;
caxis([-1 1].*max(abs(c)));
view(-37.5,10)
if save_plots, print('-r300', '-dpng', 'ridgelet_scal_Jmin03.png'); end

for j = J_min:J

   figure
   ssht_plot_sphere(plot_component(ridgelet_wav{j-J_min+1}), L, ...
      'Lighting', sphere_plot_lighting, ...
      'Type', sphere_plot_type, 'ParametricScale', sphere_plot_scale, ...
      'ParametricMin', true, ...
      'Method', sampling_method, ...
      'ColourBar', true);
   c = caxis;
   caxis([-1 1].*max(abs(c)));
   view(-37.5,10)
   
   plot_filename = sprintf('plots/ridgelet_wav_j%2.2d.png', j);
   if save_plots, print('-r300', '-dpng', plot_filename); end
   
end

