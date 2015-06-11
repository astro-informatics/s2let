function s2let_plot_cur_tiling(B, L, J_min, varargin)
% plot_cur_tiling
%  - Plot the tiling of scaling function and curvelets in harmonic space.
%  - then plot the scaling function and curvelets in real space. 
%
% Default usage :
%
%   s2let_plot_cur_tiling(B, L, J_min, <options>)
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
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

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
N = L; 


<<<<<<< HEAD
% ---------------
% Tile curvelets:
% ---------------
[cur_lm scal_l] = s2let_curvelet_tiling(args.B, args.L, args.J_min, ...
                                        'Spin', args.Spin, 'SpinLowered', args.SpinLowered,...
                                        'SpinLoweredFrom',args.SpinLoweredFrom);
% 
el_min = max(abs(args.Spin), abs(args.SpinLoweredFrom));
kappa0_cur = zeros(1,L);
for el = el_min:L-1
 kappa0_cur(1,el+1) = scal_l(el^2+el+1,1)/sqrt((2*el+1)/(4.0*pi)) ;
=======
% Tile curvelet and the scaling functions: 
[cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min);

J = s2let_jmax(L, B);
xi = 0:0.01:L-1;
x = 0:L-1;

%
% Plot the scaling function: 
%
kappa0 = zeros(1,L);
for l = 0:L-1
 kappa0(1,l+1) = scal_l(l^2+l+1,1);
>>>>>>> parent of cc46a19... 5th commit: add transform_mw.m and demo 8 and 9
end

kappa = zeros(J+1,L);
for j = J_min:J
<<<<<<< HEAD
 for el = el_min:L-1 % 0:L-1
=======
 for l = 0:L-1
>>>>>>> parent of cc46a19... 5th commit: add transform_mw.m and demo 8 and 9
  % ind = l^2 +l + m + 1 ; now consider m =  el; 
  kappa(j+1,l+1) = cur_lm{j-J_min+1}(1,l^2+l+l+1);
 end
end 


<<<<<<< HEAD
% Set for the output figures: 
pltroot = '../../../figs/' ;
configstr = ['Spin',int2str(args.Spin),...
             '_N',int2str(N),'_L',int2str(L),'_B',int2str(B),...
             '_Jmin',int2str(J_min)];


% 
xi =0:0.01:L-1;
x = 0:L-1;
% ------------
% Plot the tiling of the scaling function: 
% ------------
=======
>>>>>>> parent of cc46a19... 5th commit: add transform_mw.m and demo 8 and 9
figure('Position',[100 100 900 450])
  %semilogx(0:L-1, kappa0, 'k', 'LineWidth', 2);
yi = interp1(x, kappa0, xi,'pchip');
semilogx(xi, yi, 'k', 'LineWidth', 2);
  %h = text(2, 1.07, 'k0', 'Color', [0 0 0]);
hold on;
for j = J_min:J
  colour = rand(1,3)*0.9;
  %plot(0:L-1, kappa(j+1,:), 'LineWidth', 2, 'Color', colour);
  yi = interp1(x,kappa(j+1,:),xi,'pchip');
  plot(xi, yi, 'LineWidth', 2, 'Color', colour);
  %h = text(B.^j, 1.07, strcat('j',num2str(j+1)), 'Color', colour);
end
%title('Harmonic tiling');
%xlabel('l');
axis([0 L -0.05 1.15]);
set(gca,'XTick',2.^[0:(J+2)]);


%
% Plot the scaling function in real space: 
%
nx = 2;
ny = 3;
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L);
figure('Position',[100 100 900 200]) 
h = subplot(nx, ny, 1);
f = ssht_inverse(scal_l, L, 'Reality', true);
plot(thetas, f(:,1), '-k', 'LineWidth', 2)
mx = 1.1*max(f(:,1));
axis([0 3. -mx/8 mx ]) 


%
% Plot the curvelet kernels in real space:  
%
Jmax = J;
for j = J_min:Jmax
   h = subplot(nx, ny, j-J_min+2);
   hold on
   flm = zeros(L^2,1);
    for l = 0:L-1
        flm(l^2+l+1,1) = kappa(j+1,l+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), '-k', 'LineWidth', 2) 
   mx = 1.1*max(f(:,1));
   axis([0 3. -mx/1.5 mx ])
end 


end
