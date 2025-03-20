% s2let_demo5
% Compute and plot the harmonic tiling and the wavelet kernels, and output png figures.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

B = 2;
zoomfactor = 1.4;
J_min = 2;
L = 64;
Spin = 0;
N = 0;
J = s2let_jmax(L, B);
plot_caxis_scale = 20
type = 'colour';
lighting = false;

ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 4;
nx = 3;

maxfigs = nx*ny;
pltroot = '../../../figs'
configstr = ['Spin',int2str(Spin),'_N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]
[psi_lm phi_l] = s2let_wavelet_tiling(B, L, N, Spin, J_min);

figure('Position',[100 100 300 300])
h=subplot(1, 1, 1);
flm = zeros(L^2,1);
for l = 0:L-1
    flm(l^2+l+1,1) = phi_l(l+1);
end
f = ssht_inverse(flm, L, 'Reality', true);
ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
title(h,'Scaling fct')
locate = get(h,'title');
pos = get(locate,'position');
pos(1,2) = pos(1,2)+0.7;
pos(1,1) = pos(1,1)-0.7;
set(locate,'pos',pos);
zoom(1.2)%
%v = caxis;
%temp = max(abs(v));
%caxis([-temp temp]*plot_caxis_scale)


%colormap(jet)
%fname = [pltroot,'/s2let_demo5_', configstr, '_scal_jet.png']
%print('-r200', '-dpng', fname)
%colormap(hot)
%fname = [pltroot,'/s2let_demo5_', configstr, '_scal_hot.png']
%print('-r200', '-dpng', fname)

figure('Position',[100 100 1200 600])
ind = 0;
for j = J_min:J
   flm = psi_lm(:,j+1);

   if Spin == 0
       f = ssht_inverse(flm, L, 'Reality', true);
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
           ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
           title(h, ['Wavelet j = ',int2str(j-J_min+1)])
           locate = get(h,'title');
           pos = get(locate,'position');
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           %v = caxis;
           %temp = max(abs(v));
           %caxis([-temp temp]*plot_caxis_scale)
           %zoom(zoomfactor)
       end
   end
   if Spin > 0
       f = ssht_inverse(flm, L, 'spin', Spin);
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
           ssht_plot_sphere(real(f), L, 'Type', type, 'Lighting', lighting);
           title(h, ['Wavelet j = ',int2str(j-J_min+1), ', real part'])
           locate = get(h,'title');
           pos = get(locate,'position');
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           %v = caxis;
           %temp = max(abs(v));
           %caxis([-temp temp]*plot_caxis_scale)
           %zoom(zoomfactor)
       end
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
           ssht_plot_sphere(imag(f), L, 'Type', type, 'Lighting', lighting);
           title(h, ['Wavelet j = ',int2str(j-J_min+1), ', imag part'])
           locate = get(h,'title');
           pos = get(locate,'position');
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           %v = caxis;
           %temp = max(abs(v));
           %caxis([-temp temp]*plot_caxis_scale)
           %zoom(zoomfactor)
       end
       
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
           ssht_plot_sphere(abs(f), L, 'Type', type, 'Lighting', lighting);
           title(h, ['Wavelet j = ',int2str(j-J_min+1), ', abs part'])
           locate = get(h,'title');
           pos = get(locate,'position'); 
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           %v = caxis;
           %temp = max(abs(v));
           %caxis([-temp temp]*plot_caxis_scale)
           %zoom(zoomfactor)
       end
       
   end
end


colormap(jet)
fname = [pltroot,'/s2let_demo5_', configstr, '_wav_jet.png']
print('-r200', '-dpng', fname)
colormap(hot)
fname = [pltroot,'/s2let_demo5_', configstr, '_wav_hot.png']
print('-r200', '-dpng', fname)

