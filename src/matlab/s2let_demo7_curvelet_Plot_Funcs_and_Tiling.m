% s2let_demo7 - plotting curvelets 
%
% 1) Plot the scaling function and curvelet functions in real space (s2let_plot_curvelet_tiling).
% 2) Plot the curvelets on the Sphere (s2let_plot_curvelet_on_sphere).
% 3) Plot the curvelets parametrically (s2let_plot_curvelet_parametric).
% wherein  curvelets and the scaling functions are generated via
% the matlab function "s2let_curvelet_tiling(B, L, J_min, <options>)".
%
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

clear all;

% ---------------
% Define curvelet parameters: 
% ---------------
Spin = 0;            % Spin value of curvelet 
B = 2;               % B = 2 for dyadic sampling
L = 128;             % Angular band-limit
J_min = 2;           % Minimum scale probed by curvelets
J =s2let_jmax(L, B); % Maximum scale probed by curvelets =ceil(log L/ log B);


% ---------------
% Plot the tiling of scaling function and curvelets in harmonic space:
% And plot the scaling function and curvelet functions in real space:
% ---------------
disp(' - Plot the tiling of curvelets');
s2let_plot_curvelet_tiling(B, L, J_min,...
                           'Spin', Spin);

% Define Euler angles for rotating the curvelets functions 
% such that the curvelets centred on the North pole of the sphere: 
alpha =  0 ;
beta = acos(-Spin/L); 
gamma = 0  ;  

% ---------------
% Plot the curvelets on the sphere:
% ---------------     
disp(' - Plot the curvelets on the sphere');
s2let_plot_curvelet_on_sphere(alpha, beta, gamma, B, L, J_min, ...
                              'Spin', Spin);

%{
disp(' - Plot the curvelet parametrically');
s2let_plot_curvelet_parametric(alpha, beta, gamma,...
                               B, L, J_min, ...
                               'Spin', Spin);
%}
