% s2let_demo7_curvelet_Plot_Funcs_and_Tiling
%
% 1) Plot the scaling function and curvelet functions in real space (s2let_plot_curvelet_tiling).
% 2) Plot the curvelets on the Sphere (s2let_plot_curvelet_on_sphere).
% wherein  curvelets and the scaling functions are generated via
% the matlab function "s2let_curvelet_tiling(B, L, J_min, <options>)".
%
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2015  Boris Leistedt, Martin Büttner,
%                     Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------


clear all;

% ---------------
% Define curvelet parameters: 
% ---------------
Spin = 0;            % Spin value of curvelet 
B = 2;               % B = 2 for dyadic sampling
L = 64;              % Angular band-limit
N= L;                % for curvelet, the azimuthal band-limit N = L
J_min = 2;           % Minimum scale probed by curvelets
J =s2let_jmax(L, B); % Maximum scale probed by curvelets =ceil(log L/ log B);


% ---------------
% Plot the tiling of scaling function and curvelets in harmonic space:
% And plot the scaling function and curvelet functions in real space:
% ---------------
disp(' - Plot the tiling of curvelets');
s2let_plot_curvelet_tiling(B, L, J_min,...
                           'Spin', Spin);

% ---------------
% Plot the curvelets on the sphere:
% ---------------     
% Define Euler angles for rotating the curvelets functions 
% such that the curvelets centred on the North pole of the sphere (instead of along -x direction):
alpha =  0 ;
beta = pi/2 ;
gamma = 0 ;
disp(' - Plot the curvelets on the sphere');
s2let_plot_curvelet_on_sphere(alpha, beta, gamma,...
                              B, L, J_min, ...
                              'Spin', Spin);

