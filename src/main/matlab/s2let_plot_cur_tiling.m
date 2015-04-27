function s2let_plot_cur_tiling(B, L, J_min)

% OR: s2let_plot_cur_tiling(B, L, N, Spin, J_min, varargin)

% plot_cur_tiling - Plot tiling in harmonic space.
% -- Curvelets on the sphere.
%
% Default usage :
%
%   plot_cur_tiling(B, L, J_min)
%
% B is the wavelet parameter,
% L is the angular band-limit,
% J_min the first wavelet to be used.
%
% N is the azimuthal/directional band-limit,
% Spin is the spin number,
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'SpinLowered' = { true  [Apply normalisation factors for spin-lowered
                            %                           wavelets and scaling function.],
    %                    false [Apply the usual normalisation factors such
                                %                           that the wavelets fulfil the admissibility
                                %                           condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
                        %                       option indicates which spin number the wavelets
                        %                       should be lowered from (default = 0)]
%
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details


[kappa kappa0] = s2let_transform_axisym_tiling(B, L, J_min);
%[psi_lm phi_l]  = s2let_curvelet_tiling(B, L, N, Spin, J_min, varargin);

J = s2let_jmax(L, B);
xi = 0:0.01:L-1;
x = 0:L-1;

figure('Position',[100 100 900 450])
%semilogx(0:L-1, kappa0, 'k', 'LineWidth', 2);
yi = interp1(x,kappa0,xi,'pchip');
semilogx(xi, yi, 'k', 'LineWidth', 2);
%h = text(2, 1.07, 'k0', 'Color', [0 0 0]);
hold on;
for j = J_min:J
  colour = rand(1,3)*0.9;
  %plot(0:L-1, kappa(j+1,:), 'LineWidth', 2, 'Color', colour);
    yi = interp1(x,kappa(j+1,:),xi,'pchip');
    plot(xi, yi, 'LineWidth', 2, 'Color', colour);
  %h = text(B.^j, 1.07, strcat('j',num2str(j+1)), 'Color', colour);
%
% Directional Wavelet - load('kappas_s2dw'):
% yi = interp1(x, kappa_s2dw(j+1,:), xi,'pchip');
% plot(xi, yi, '-k', 'LineWidth', 2)%, 'Color', colour);
%
end
%title('Harmonic tiling');
%xlabel('el');
axis([0 L -0.05 1.15]);
set(gca,'XTick',2.^[0:(J+2)]);

end
