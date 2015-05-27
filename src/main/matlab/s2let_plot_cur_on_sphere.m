function s2let_plot_cur_on_sphere(alpha, beta, gamma, cur_lm, scal_l, L, J_min, varargin)
%
% s2let_plot_cur_on_sphere -
% Plot curvelet coefficients on multiple spheres.
%
% This function
% i) compute the j-th curvelet, rotated by rho=(alpha, beta, gamma) in
% harmonic space and reconstruct it on the sphere.
% ii) generates one plot of the scaling function contribution and
% a grid of plots for each orientation of each scale of the
% curvelet contributions.
%
% Default usage :
%
%   s2let_plot_cur_on_sphere(alpha, beta, gamma, cur_lm, scal_l, L, J_min, <options>)
%
% cur is cell array with all the curvelet coefficients.
% its first index is the wavelet scale j, the second
% index is the orientation g, and each element is a
% function on the sphere in MW sampling.
% scal is the corresponding scaling function contribution
% (i.e. just a single function on the sphere).
% B is the wavelet parameter.
% L is the angular band-limit.
% N is the orientational band-limit.
% J_min is the first curvelet scale in cur.
%
% j is the order of the curvelet under consideration (depends on B)
% rho=(alpha, beta, gamma) is the rotation in SO(3) by which to rotate
% the curvelet
% L if harmonic band-limit for the reconstruction on the sphere
% psi_j is the reconstructed curvelet on the sphere, at resolution L
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'B'               = { Dilation factor; B > 1 (default = 2) }
%  'N'               = { Azimuthal band-limit; N > 0 (default = L) }
%  'Spin'            = { Spin number; Spin >= 0 (default = 0) }
%  'Upsample'      = { false   [multiresolution algorithm (default)],
%                      true  [full resolution wavelets] },
%  'Function'        = { 'real' [plot the real part of the input functions (default)],
%                        'imag' [plot the imaginary part of the input functions],
%                        'abs'  [plot the absolute value of the input functions] }
%
% S2LET package to perform curvelet transform on the Sphere.
% Copyright (C) 2012-2014  Boris Leistedt, Martin Büttner & Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('alpha', @isnumeric);
p.addRequired('beta', @isnumeric);
p.addRequired('gamma', @isnumeric);
p.addRequired('cur_lm', @iscell);
p.addRequired('scal_l', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('N', -1, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Function', 'real', @ischar)
p.parse(alpha, beta, gamma, cur_lm, scal_l, L, J_min, varargin{:});


args = p.Results;

if args.N == -1
    args.N = L;
end

B = args.B;
N = args.N;
Spin = args.Spin;
alpha = args.alpha;
beta = args.beta;
gamma = args.gamma ;
J_min = args.J_min;
J = s2let_jmax(L, B);

% Precompute Wigner small-d functions
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end

% Define plotting parameters
zoomfactor = 1.4;
plot_caxis_scale = 2
type = 'colour';
lighting = true;
nx = 3;
ny = 4;
maxfigs = nx*ny;
pltroot = '../../../figs/' ;
configstr = ['Spin',int2str(Spin),'_N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)];

%
% curvelets:
%
figure('Position',[100 100 1200 600])
ind=0;
for j = J_min:J,
%% Rotate the curvelets coefficients
flm_cur_rot = ssht_rotate_flm(cur_lm{j-J_min+1}, d, alpha, gamma);
  if Spin == 0
   % Compute the function (rotated):
   f_cur_rot = ssht_inverse(flm_cur_rot, L, 'Reality', true);
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

   if Spin > 0
    f_cur_rot = ssht_inverse(flm_cur_rot, L, 'spin', Spin);
    ind = ind + 1;
     if ind <= maxfigs
      h = subplot(ny, nx, ind);
      ssht_plot_sphere(real(f_cur_rot), L, 'Type', type, 'Lighting', lighting);
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
     ssht_plot_sphere(abs(f_cur_rot), L, 'Type', type, 'Lighting', lighting);
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
fname = [pltroot,'s2let_plotfn_', configstr, '_cur_jet.png'];
print('-r200', '-dpng', fname);


%
%  Scaling functions
%
figure('Position',[100 100 300 300])
h=subplot(1, 1, 1);
f = ssht_inverse(scal_l, L, 'Reality', true);
ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
% f_scal = ssht_inverse(scal_l, L, 'Reality', true);
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
fname = [pltroot,'s2let_plotfn_', configstr, '_scal_jet.png'] ;
print('-r200', '-dpng', fname);

