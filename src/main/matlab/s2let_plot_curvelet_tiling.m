function s2let_plot_curvelet_tiling(B, L, J_min, varargin)
% s2let_plot_curvelet_tiling
%  - Plot the tiling of scaling functions and curvelets in harmonic space.
%  - then plot the scaling functions and curvelets in real space. 
%
% Default usage :
%
%   s2let_plot_curvelet_tiling(B, L, J_min, <options>)
%
% B is the curvelet dilation parameter,
% L is the angular band-limit,
% J_min the first wavelet to be used.
%
% % Valid options include:
%
%  'Spin'        = { Spin; (default=0) }
%  'SpinLowered' = { true  [Apply normalisation factors for spin-lowered
%                           wavelets and scaling function.],
%                    false [Apply the usual normalisation factors such
%                           that the wavelets fulfil the admissibility
%                           condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the wavelets
%                       should be lowered from (default = 0)]
%
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

% Parse arguments.
p = inputParser;
p.addRequired('B', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(B, L, J_min, varargin{:});

args = p.Results;

B = args.B;
L = args.L;
J_min = args.J_min;
Spin = args.Spin;
J = s2let_jmax(L, B);

% ---------------
% Tile curvelets:
% ---------------
[cur_lm scal_l] = s2let_curvelet_tiling(args.B, args.L, args.J_min, ...
                                        'Spin', args.Spin, 'SpinLowered', args.SpinLowered,...
                                        'SpinLoweredFrom',args.SpinLoweredFrom);
                                    
% Normalise and reshape the scaling functions: 
el_min = max(abs(args.Spin), abs(args.SpinLoweredFrom));
kappa0_cur = zeros(1,L);
for el = el_min:L-1
 kappa0_cur(1,el+1) = scal_l(el^2+el+1,1)/sqrt((2*el+1)/(4.0*pi)) ;
end 

% Normalise and reshape the curvelet functions: 
kappa_cur = zeros(J+1,L);
for j = J_min:J
 for el= el_min:L-1
  % ind = l^2 +l + m + 1 ; now consider m =  el; 
  kappa_cur(j+1,el+1) = cur_lm{j-J_min+1}(1,el^2+el+el+1)/ ...
                        (sqrt(1./2.)* sqrt((2*el+1)/(8.0*pi*pi))) ;
 end
end 


% Set for the output figures: 
pltroot = '../../../figs/' ;
configstr = ['Spin',int2str(args.Spin),...
             '_L',int2str(L),'_B',int2str(B),...
             '_Jmin',int2str(J_min)];


% 
xi =0:0.01:L-1;
x = 0:L-1;
% ------------
% Plot the tiling of the scaling function: 
% ------------
figure('Position',[100 100 900 450])
yi = interp1(x, kappa0_cur, xi,'pchip');
semilogx(xi, yi, 'k', 'LineWidth', 2);
  %h = text(2, 1.07, 'k0', 'Color', [0 0 0]);
hold on;
% ------------
% Plot the tiling of the curvelet kernels : 
% ------------
for j = J_min:J
  colour = rand(1,3)*0.9;
  yi = interp1(x,kappa_cur(j+1,:),xi,'pchip');
  plot(xi, yi, 'LineWidth', 2, 'Color', colour);
  %h = text(B.^j, 1.07, strcat('j',num2str(j+1)), 'Color', colour);
end

%title('Harmonic tiling');
%xlabel('l');
axis([0 L -0.05 1.15]);
set(gca,'XTick',2.^[0:(J+2)]);
fname = [pltroot,'s2let_plot_cur_tiling_', configstr, '.png']
print('-r200', '-dpng', fname);



% ----------------------------
% Plot the scaling function in real space: 
% ----------------------------
nx = 3;
ny = 3;
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L);
figure('Position',[100 100 900 200]) 
h = subplot(nx, ny, 1);
f = ssht_inverse(scal_l, L, 'Reality', true);
plot(thetas, f(:,1), '-k', 'LineWidth', 2)
mx = 1.1*max(f(:,1));
axis([0 3. -mx/8 mx ]) 

% ----------------------------
% Plot the curvelet kernels in real space:  
% ----------------------------
Jmax = J;
for j = J_min:Jmax
   h = subplot(nx, ny, j-J_min+2);
   hold on
   flm = zeros(L^2,1);
    for el = el_min:L-1
        flm(el^2+el+1,1) = kappa_cur(j+1,el+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), '-k', 'LineWidth', 2) 
   mx = 1.1*max(f(:,1));
   axis([0 3. -mx/1.5 mx ])
end 
fname = [pltroot,'s2let_plot_cur_tiling_fn_real', configstr, '.png']
print('-r200', '-dpng', fname);

end

