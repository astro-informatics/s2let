% s2let_plot_wavs

load('kappas_spline');
load('kappas_s2dw');
load('kappas_need');

B = 3;
J_min = 2;
L = 128;
J = s2let_jmax(L, B);
Jmax = 3;

ns = ceil(sqrt(2+J-J_min+1)) ;
nx = 1;
ny = 3;



J = s2let_jmax(L, B);
xi = 0:0.01:L-1;
x = 0:L-1;


figure('Position',[100 100 900 450])
yi = interp1(x,kappa0_spline,xi,'pchip');
semilogx(xi, yi, '-.r', 'LineWidth', 2);
hold on;
yi = interp1(x,kappa0_s2dw,xi,'pchip');
plot(xi, yi, '-k', 'LineWidth', 2);
yi = interp1(x,kappa0_need,xi,'pchip');
plot(xi, yi, '--b', 'LineWidth', 2);
for j = J_min:J  
  colour = rand(1,3)*0.9;
  yi = interp1(x, kappa_spline(j+1,:), xi,'pchip');
  plot(xi, yi, '-.r', 'LineWidth', 2)%, 'Color', colour);
  yi = interp1(x, kappa_s2dw(j+1,:), xi,'pchip');
  plot(xi, yi, '-k', 'LineWidth', 2)%, 'Color', colour);
  yi = interp1(x, kappa_need(j+1,:), xi,'pchip');
  plot(xi, yi, '--b', 'LineWidth', 2)%, 'Color', colour);
end
axis([1 L -0.05 1.15]);
set(gca,'XTick',2.^[0:(J+2)]);
hleg1 = legend('B-Spline', 'SD', 'Needlet');
set(hleg1, 'Position', [.15,.25,.1,.2]);


[thetas, phis, n, ntheta, nphi] = ssht_sampling(L);
figure('Position',[100 100 900 200]) 

h = subplot(nx, ny, 1);
hold on
flm = zeros(L^2,1);
for l = 0:L-1
    flm(l^2+l+1,1) = kappa0_spline(l+1);
end     
f = ssht_inverse(flm, L, 'Reality', true);
plot(thetas, f(:,1), '-.r', 'LineWidth', 2)
mx = 1.1*max(f(:,1));
axis([0 2. -mx/8 mx ])
flm = zeros(L^2,1);
for l = 0:L-1
    flm(l^2+l+1,1) = kappa0_s2dw(l+1);
end     
f = ssht_inverse(flm, L, 'Reality', true);
plot(thetas, f(:,1), '-k', 'LineWidth', 2)
   
flm = zeros(L^2,1);
for l = 0:L-1
    flm(l^2+l+1,1) = kappa0_need(l+1);
end     
f = ssht_inverse(flm, L, 'Reality', true);
plot(thetas, f(:,1), '--b', 'LineWidth', 2)
   
   
for j = J_min:Jmax
   h = subplot(nx, ny, j-J_min+2);
   hold on
   flm = zeros(L^2,1);
    for l = 0:L-1
        flm(l^2+l+1,1) = kappa_spline(j+1,l+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), '-.r', 'LineWidth', 2) 
    mx = 1.1*max(f(:,1));
    axis([0 2. -mx/7 mx ])
   flm = zeros(L^2,1);
    for l = 0:L-1
        flm(l^2+l+1,1) = kappa_s2dw(j+1,l+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), '-k', 'LineWidth', 2) 
   flm = zeros(L^2,1);
    for l = 0:L-1
        flm(l^2+l+1,1) = kappa_need(j+1,l+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   plot(thetas, f(:,1), '--b', 'LineWidth', 2) 
end
