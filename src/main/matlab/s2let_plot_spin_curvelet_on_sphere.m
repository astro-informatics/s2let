function s2let_plot_spin_curvelet_on_sphere(B, L, J_min, varargin)
%
% s2let_plot_cur_on_sphere -
% Plot curvelet coefficients on multiple spheres.
%
% This function
% i) compute the j-th curvelet, which has been rotated by the Euler's angles (alpha =pi, beta =pi/2, gamma=2) in
% harmonic space (in s2let_spin_curvelet_tiling.m), 
% and reconstruct it on the sphere.
% ii) generates one plot of the scaling function contribution and
% a grid of plots for each orientation of each scale of the
% curvelet contributions.
%
% Default usage :
%
%   s2let_plot_curvelet_on_sphere(B, L, J_min, <options>)
%
% B is the wavelet dilation factor 
% L is the angular band-limit.
% J_min is the first curvelet scale in cur.
%
% Option :
%  'N'               = { Azimuthal band-limit; N > 0 , and should =L for curvelets; (default = L) }
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
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%
%
% j is the order of the curvelet under consideration (depends on B)
% L if harmonic band-limit for the reconstruction on the sphere
% psi_j is the reconstructed curvelet on the sphere, at resolution L
%
% ---------------------------------------------------------
% S2LET package to perform wavelet transforms on the sphere.
% Copyright (C) 2012-2014 Boris Leistedt and Jason McEwen
% See LICENSE.txt for license details
%
% Modified S2LET package to perform curvelet transforms on the Sphere.
% ---------------------------------------------------------

% Parse arguments.
p = inputParser;
p.addRequired('B', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('Sampling', 'MW', @ischar);
p.parse(B, L, J_min, varargin{:});
args = p.Results;

B = args.B ; 
L = args. L; 
N = L;
J_min = args.J_min; 
J = s2let_jmax(L, B);

[cur_lm scal_l] = s2let_spin_curvelet_tiling(B, L, J_min, ...
                                             'Spin', args.Spin, ...
                                             'SpinLowered', args.SpinLowered,...
                                             'SpinLoweredFrom', args.SpinLoweredFrom);
                                    
% Define plotting parameters
zoomfactor = 1.4;
plot_caxis_scale = 2;
type = 'colour';
lighting = true;
nx = 3;
ny = 4;
maxfigs = nx*ny;
pltroot = '../../../figs/' ;
configstr = ['Spin',int2str(args.Spin),...
             '_N',int2str(N),'_L',int2str(L),'_B',int2str(B),...
             '_Jmin',int2str(J_min)];

%
% curvelets:
%
figure('Position',[100 100 1200 600])
ind=0;
for j = J_min:J,
%% Rotate the curvelets coefficients
  if args.Spin == 0
   % Compute the curvelet functsions (from the roatated curvelet harmonic coefficients)
   % where size(f_cur_rot)= [L, 2*L-1]
   f_cur_rot = ssht_inverse(cur_lm{j-J_min+1}, L, ...
                           'Method', args.Sampling, ...
                           'Spin', args.Spin, 'Reality', true);
   ind = ind + 1;
    if ind <= maxfigs
     h = subplot(ny, nx, ind);
     %% Plot the rotated function on the sphere
     ssht_plot_sphere(f_cur_rot, L, 'Type', type, 'Lighting', lighting);
     title(h, ['Curvelet j = ',int2str(j-J_min+1)])
     locate = get(h,'title');
     pos = get(locate,'position');
     pos(1,2) = pos(1,2)+0.7;
     pos(1,1) = pos(1,1)-0.7;
     set(locate,'pos',pos);
     v = caxis;
     temp = max(abs(v));
     caxis([-temp temp]*plot_caxis_scale);
     zoom(zoomfactor)
    end
   end

   if args.Spin > 0
    % Compute the curvelet functsions (from the roatated curvelet harmonic coefficients cur_lm{j-J_min+1} whose size is (L^2,1)) 
    % where size(f_cur_rot) = [L, 2L-1]  
    f_cur_rot = ssht_inverse(cur_lm{j-J_min+1}, L, ...
                            'Method', args.Sampling,...
                            'Spin', args.Spin);            
    ind = ind + 1;
     if ind <= maxfigs
      h = subplot(ny, nx, ind);
      ssht_plot_sphere(real(f_cur_rot), L,  'Type', type,'Lighting', lighting);
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
     end

    ind = ind + 1;
    if ind <= maxfigs
     h = subplot(ny, nx, ind);
     ssht_plot_sphere(imag(f_cur_rot), L, 'Type', type, 'Lighting', lighting);
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
    end

    ind = ind + 1;
    if ind <= maxfigs
     h = subplot(ny, nx, ind);
     ssht_plot_sphere(abs(f_cur_rot), L, 'Type', type,'Lighting', lighting);
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
    end
   end
end
% output as png file
colormap(jet)
fname = [pltroot,'s2let_plotfn_', configstr, '_cur_jet.png']
print('-r200', '-dpng', fname);

%
%  Scaling functions
%
figure('Position',[100 100 300 300])
h=subplot(1, 1, 1);
f = ssht_inverse(scal_l, L, 'Reality', true);
ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
%
title(h,'Scaling function')
locate = get(h,'title');
pos = get(locate,'position');
pos(1,2) = pos(1,2)+0.7;
pos(1,1) = pos(1,1)-0.7;
set(locate,'pos',pos);
zoom(1.2)
v = caxis;
temp = max(abs(v));
caxis([-temp temp]*plot_caxis_scale)
% output as png file
colormap(jet)
fname = [pltroot,'s2let_plotfn_', configstr, '_scal_jet.png'] 
print('-r200', '-dpng', fname);

