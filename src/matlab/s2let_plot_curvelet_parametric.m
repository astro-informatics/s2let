function s2let_plot_curvelet_parametric(alpha, beta, gamma, B, L, J_min, varargin)
% s2let_plot_curvelet_parametric
% - Plot curvelet coefficients parametrically on multiple spheres.
%
% This Matlab function
% i)  compute the j-th curvelet, rotated by rho=(alpha, beta, gamma) in
%     harmonic space and reconstruct it on the sphere.
% ii) plot the scaling function contribution and curvelet coefficients parametrically
%
% Default usage :
%
%   s2let_plot_curvelet_parametric(alpha, beta, gamma, B, L, J_min, <options>)
%
% (alpha, beta, gamma) are the Euler's angles by which we rotate the curvelet
% B is the wavelet dilation factor 
% L is the angular band-limit.
% J_min is the first curvelet scale to be probed.
%
% Option :
%  'Spin'            = { Spin number; Spin >= 0 (default = 0) }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                              wavelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                              that the wavelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the wavelets
%                       should be lowered from (default = 0)]
%  'Reality'         = { false        [do not assume corresponding signal f real (default)],
%                        true         [assume f real (improves performance)] }
%
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

% Parse arguments.
p = inputParser;
p.addRequired('alpha', @isnumeric);
p.addRequired('beta', @isnumeric);
p.addRequired('gamma', @isnumeric);
p.addRequired('B', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('Sampling', 'MW', @ischar);
p.parse(alpha, beta, gamma, B, L, J_min, varargin{:});
args = p.Results;

B = args.B ; 
L = args. L;
J_min = args.J_min; 
J = s2let_jmax(L, B);

original_spin = 0 ;  % if we don't use spin-lowered wavelets (default). 
if (args.SpinLowered ~= 0) % For spin-lowered curvelet: 
    original_spin= args.SpinLoweredFrom; 
end 
el_min = max(abs(args.Spin), abs(original_spin));

% Precompute Wigner small-d functions for rotation 
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, args.beta);  %el_min, beta);
for el = 1:L-1  %el_min:L-1  %
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, args.beta);
end

% Define plotting parameters
zoomfactor = 0.5;
plot_caxis_scale = 2;
sampling_method = args.Sampling;
sphere_plot_type = 'parametric';
sphere_plot_lighting = true;
sphere_plot_scale = 1;

plot_component = @(x) real(x);

nx = 3;
ny = 5;
maxfigs = nx*ny;
pltroot = '../../../figs/' ;
configstr = ['Spin',int2str(args.Spin),...
             '_L',int2str(L),'_B',int2str(B),...
             '_Jmin',int2str(J_min)];

% ------------ Tile curvelts and the scaling functions --------------- 
[cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min, ...
                                        'Spin', args.Spin, ...
                                        'SpinLowered', args.SpinLowered,...
                                        'SpinLoweredFrom', args.SpinLoweredFrom);          

% -------------
% Plot curvelets:
% -------------
% For Spin-0 
if args.Spin == 0
    for j = J_min:J,
%% Rotate the curvelets coefficients
        flm_cur_rot = ssht_rotate_flm(cur_lm{j-J_min+1}(:), d, args.alpha, args.gamma);
        % Compute the function (rotated):
        f_cur_rot = ssht_inverse(flm_cur_rot, L, 'Method', args.Sampling, 'Spin', args.Spin, 'Reality', true);
%% Plot the rotated function on the sphere
         figure; 
         ssht_plot_sphere(plot_component(f_cur_rot), L, ...
                          'Lighting', sphere_plot_lighting, ...
                          'Type', sphere_plot_type, ...
                          'ParametricScale', sphere_plot_scale, ...
                          'ParametricMin', true, ...
                          'Method', sampling_method, ...
                          'ColourBar', true);
          title(['Curvelet j = ',int2str(j-J_min+1)])
          v = caxis;
          temp = max(abs(v));
          caxis([-temp temp]*plot_caxis_scale);
          zoom(zoomfactor)
          view(90,20)
  end
end

  
if (args.Spin ~= 0)
    ind =0 ; 
    figure; 
    for j = J_min:J,
%% Rotate the curvelets coefficients
        flm_cur_rot = ssht_rotate_flm(cur_lm{j-J_min+1}(:), d, args.alpha, args.gamma);
        f_cur_rot = ssht_inverse(flm_cur_rot, L, 'Method', args.Sampling,...
                                 'Spin', args.Spin,'Reality',  args.Reality);
        ind = ind + 1;
        h = subplot(ny, nx, ind);
        if ind <= maxfigs
           h = subplot(ny, nx, ind);
           ssht_plot_sphere(plot_component(f_cur_rot), L, ...
                            'Lighting', sphere_plot_lighting, ...
                            'Type', sphere_plot_type, ...
                            'ParametricScale', sphere_plot_scale, ...
                            'ParametricMin', true, ...
                            'Method', sampling_method, ...
                            'ColourBar', true);
     
            title(h, ['Spin Curvelet j = ',int2str(j-J_min+1), ', real part'])
            locate = get(h,'title');
            pos = get(locate,'position');
            pos(1,2) = pos(1,2)+0.7;
            pos(1,1) = pos(1,1)-0.7;
            set(locate,'pos',pos);
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp]*plot_caxis_scale)
            zoom(zoomfactor)
            view(70,26)
        end
   
% Plot imaginary part: 
        plot_component = @(x) imag(x);
        ind = ind + 1;
        if ind <= maxfigs
            h = subplot(ny, nx, ind);
            ssht_plot_sphere(plot_component(f_cur_rot), L, ...
                             'Lighting', sphere_plot_lighting, ...
                             'Type', sphere_plot_type, ...
                             'ParametricScale', sphere_plot_scale, ...
                             'ParametricMin', true, ...
                             'Method', sampling_method, ...
                             'ColourBar', true);
            title(h, ['Spin Curvelet j = ',int2str(j-J_min+1), ', imag part'])
            locate = get(h,'title');
             pos = get(locate,'position');
             pos(1,2) = pos(1,2)+0.7;
             pos(1,1) = pos(1,1)-0.7;
             set(locate,'pos',pos);
             v = caxis;
             temp = max(abs(v));
             caxis([-temp temp]*plot_caxis_scale)
             zoom(zoomfactor)
             view(70,26)
        end
 % Plot absolute part: 
             plot_component = @(x) abs(x);
             ind = ind + 1;
             if ind <= maxfigs
                 h = subplot(ny, nx, ind);
                 ssht_plot_sphere(plot_component(f_cur_rot), L, ...
                 'Lighting', sphere_plot_lighting, ...
                 'Type', sphere_plot_type, ...
                 'ParametricScale', sphere_plot_scale, ...
                 'ParametricMin', true, ...
                 'Method', sampling_method, ...
                 'ColourBar', true);
                 title(h, ['Spin Curvelet j = ',int2str(j-J_min+1), ', abs part'])
                 locate = get(h,'title');
                 pos = get(locate,'position');
                 pos(1,2) = pos(1,2)+0.7;
                 pos(1,1) = pos(1,1)-0.7;
                 set(locate,'pos',pos);
                 v = caxis;
                 temp = max(abs(v));
                 caxis([-temp temp]*plot_caxis_scale)
                 zoom(zoomfactor)
                 view(70,26)
             end
    end  
end
% output as png file
colormap(jet)
fname = [pltroot,'s2let_plotfn_', configstr, '_cur_jet.png']
print('-r200', '-dpng', fname);

% -------------
% Plot scaling functions
% -------------
plot_component = @(x) real(x);
figure('Position',[100 100 300 300])
h=subplot(1, 1, 1);
f_scal = ssht_inverse(scal_l, L, 'Reality', true);
ssht_plot_sphere(plot_component(f_scal), L, ...
   'Lighting', sphere_plot_lighting, ...
   'Type', sphere_plot_type, 'ParametricScale', sphere_plot_scale, ...
   'Method', sampling_method, ...
   'ColourBar', true);
%
title(h,'Scaling function')
zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp]*plot_caxis_scale)
view(0,-10)
% output as png file
colormap(jet)
fname = [pltroot,'s2let_plotfn_', configstr, '_scal_jet.png'] 
print('-r200', '-dpng', fname);

