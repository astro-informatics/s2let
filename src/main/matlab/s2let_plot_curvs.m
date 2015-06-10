% s2let_plot_curvs
% ------ About this demo -------- 
% 1) Plot the tilings for needlets, spline wavelet, directional wavelets and curvelets. 
% 2) Plot the functions of these different types of wavelets in real space.
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2015  
% -----------------------------------------------------------
% S2LET package to perform curvelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details


% For different types of wavlets:
load('kappas_spline');
load('kappas_s2dw');
load('kappas_need');

B = 3;
J_min = 2;
L = 128;
J = s2let_jmax(L, B);
Jmax = 3;
Spin = 0;

% For Curvelets: 
% Tile curvelet and the scaling functions: 
[cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min);
% The scaling function: 
kappa0_cur = zeros(1,L);
for el = 0:L-1
 % Normalise  
 scal_l(el^2+el+1,1)= scal_l(el^2+el+1,1)/sqrt((2*el+1)/(4.0*pi));
 kappa0_cur(1,el+1) = scal_l(el^2+el+1,1) ;
end
% The curvelet functions: 
kappa_cur = zeros(J+1,L);
for j = J_min:J
 for el = 0:L-1
  % ind = el^2 +el + m + 1 ; now consider m =  el; 
   cur_lm{j-J_min+1}(1,el^2+el+el+1)= cur_lm{j-J_min+1}(1,el^2+el+el+1) /(sqrt(1./2.)* sqrt((2*el+1)/(8.0*pi*pi)));
  kappa_cur(j+1,el+1) = cur_lm{j-J_min+1}(1,el^2+el+el+1) ;
 end
end 


% -----
% Next : plot all types of wavelets together: 
% -----
xi = 0:0.01:L-1;
x = 0:L-1;

ns = ceil(sqrt(2+J-J_min+1)) ;
nx = 1;
ny = 3;

% Set for the output figures: 
pltroot = '../../../figs/' ;
configstr = ['Spin',int2str(Spin),...
             '_L',int2str(L),'_B',int2str(B),...
             '_Jmin',int2str(J_min)];


% Scaling Function
figure('Position',[100 100 900 450])
yi = interp1(x,kappa0_spline,xi,'pchip');
semilogx(xi, yi, '-.r', 'LineWidth', 2);
hold on;
yi = interp1(x,kappa0_s2dw,xi,'pchip');
plot(xi, yi, '-k', 'LineWidth', 2);
yi = interp1(x,kappa0_need,xi,'pchip');
plot(xi, yi, '--b', 'LineWidth', 2);
yi = interp1(x, kappa0_cur, xi,'pchip');
plot(xi, yi, 'Color',[.7 .5 0], 'LineWidth', 2);
% Wavelets: 
for j = J_min:J  
  colour = rand(1,3)*0.9;
  yi = interp1(x, kappa_spline(j+1,:), xi,'pchip');
  plot(xi, yi, '-.r', 'LineWidth', 2)%, 'Color', colour);
  yi = interp1(x, kappa_s2dw(j+1,:), xi,'pchip');
  plot(xi, yi, '-k', 'LineWidth', 2)%, 'Color', colour);
  yi = interp1(x, kappa_need(j+1,:), xi,'pchip');
  plot(xi, yi, '--b', 'LineWidth', 2)%, 'Color', colour);
  yi = interp1(x, kappa_cur(j+1,:), xi,'pchip');
  plot(xi, yi, 'Color',[.7 .5 0], 'LineWidth', 2)%, 'Color', colour);
end
axis([1 L -0.05 1.15]);
set(gca,'XTick',2.^[0:(J+2)]);
hleg1 = legend('B-Spline', 'SD', 'Needlet', 'Curvelet');
set(hleg1, 'Position', [.15,.25,.1,.2]);
%
fname = [pltroot,'s2let_plot_curvs_diffTypes_tiling_', configstr, '.png']
print('-r200', '-dpng', fname);


% Plot different types of wavelets in real space: 
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L);
figure('Position',[100 100 900 200]) 

% Scaling functions:
h = subplot(nx, ny, 1);
hold on
% Spline - scaling function: 
flm = zeros(L^2,1);
for el = 0:L-1
    flm(el^2+el+1,1) = kappa0_spline(el+1);
end     
f = ssht_inverse(flm, L, 'Reality', true);
plot(thetas, f(:,1), '-.r', 'LineWidth', 2)
mx = 1.3*max(f(:,1));
axis([0 2. -mx/8 mx ])
% Directional wavelets - scaling function: 
flm = zeros(L^2,1);
for el = 0:L-1
    flm(el^2+el+1,1) = kappa0_s2dw(el+1);
end     
f = ssht_inverse(flm, L, 'Reality', true);
plot(thetas, f(:,1), '-k', 'LineWidth', 2)
% Needlets - scaling function:   
flm = zeros(L^2,1);
for el = 0:L-1
    flm(el^2+el+1,1) = kappa0_need(el+1);
end     
f = ssht_inverse(flm, L, 'Reality', true);
plot(thetas, f(:,1), '--b', 'LineWidth', 2)
% Curvelets - scaling function:   
f = ssht_inverse(scal_l, L, 'Reality', true);
plot(thetas, f(:,1), 'Color',[.7 .5 0], 'LineWidth', 2)
mx = 1.3*max(f(:,1));
axis([0 3. -mx/8 mx ]) 

% Wavelets: 
for j = J_min:Jmax
   h = subplot(nx, ny, j-J_min+2);
   hold on
   % Spline: 
   flm = zeros(L^2,1);
    for el = 0:L-1
        flm(el^2+el+1,1) = kappa_spline(j+1,el+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), '-.r', 'LineWidth', 2) 
   mx = 1.3*max(f(:,1));
   axis([0 3. -mx/8 mx ])
   % directional wavelets:
   flm = zeros(L^2,1);
    for el = 0:L-1
        flm(el^2+el+1,1) = kappa_s2dw(j+1,el+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), '-k', 'LineWidth', 2) 
   % Needlets:
   flm = zeros(L^2,1);
    for el = 0:L-1
        flm(el^2+el+1,1) = kappa_need(j+1,el+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), '--b', 'LineWidth', 2) 
   % Curvelets:
   flm = zeros(L^2,1);
    for el = 0:L-1
        flm(el^2+el+1,1) = kappa_cur(j+1,el+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), 'Color',[.7 .5 0], 'LineWidth', 2)  
   mx = 1.3*max(f(:,1));
   axis([0 3. -mx/8 mx ])
end
%
fname = [pltroot,'s2let_plot_curvs_diffTypes_fn_real', configstr, '.png']
print('-r200', '-dpng', fname);


