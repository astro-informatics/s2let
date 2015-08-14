% s2let_demo9_Spin_Curvelet_Mollweide_EarthTomography
% -----------------------------------------------------------
% Plot curvelet coefficients on multiple Mollweide projections.
% The function generates one plot of the scaling function
% contribution and a grid of plots for each orientation of
% each scale of the wavelet contributions. 
% -----------------------------------------------------------
% wav is cell array with all the wavelet coefficients.
% its first index is the wavelet scale j, the second
% index is the orientation g, and each element is a
% function on the sphere in MW sampling.
% scal is the corresponding scaling function contribution
% (i.e. just a single function on the sphere).
% B is the wavelet parameter.
% L is the angular band-limit.
% N is the orientational band-limit.
% J_min is the first wavelet scale in wav.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%  'Spin'            = { non-negative integers (default=0) }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'Upsample'        = { false        [multiresolution algorithm (default)],
%                        true       [full resolution curvelets] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               curvelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the curvelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the curvelets
%                       should be lowered from (default = 0)]
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%
%
% FOR ssht_plot_mollweide(f, L, <options>)
% where f is the sampled function and L is the harmonic band-limit.
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Method'          = { 'MW'         [McEwen & Wiaux sampling (default)],
%                        'MWSS'       [McEwen & Wiaux symmetric sampling],
%                        'DH'         [Driscoll & Healy sampling],
%                        'GL'         [Gauss-Legendre sampling] }
%  'ColourBar'       = { false        [do not add colour bar (default)],
%                        true         [add colour bar] }
%  'Mode'            = { 0            Plot amplitude, or modulus is f complex (default),
%                        1            Plot real part,
%                        2            Plot imaginaty part,
%                        3            Plot modulus and arrows for real/img angle }
%  'Spin'            = { non-negative integers (default=0) }
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2015  
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% 
% Modified S2LET package to perform Curvelets on the Sphere.
% -----------------------------------------------------------
% Remark: 
% Use default mode setting for ssht_plot_mollweide 
% (i.e. Mode 0 - ploting amplitude/ modulus of the complex signals ) 

 
clear; % all;
% close all ;

% ---------------
% Set the reality and upsample flag for signal analysis:
% ---------------
reality = false;    % true for real data;
upsample = true ;  % true for full resolution plot

% ---------------
% Load data set: 
% ---------------
%filename= '../../../data/heal_grace.fits';
filename= '../../../data/heal_galileo.fits';
%filename= '../../../data/heal_rnl.fits';
%filename= '../../../data/heal_uffizi.fits';
%filename= '../../../data/heal_stpeters.fits';

info =fitsinfo(filename); 
disp(info.Contents);
info.BinaryTable
rowend = info.BinaryTable.Rows;   

figure

L = 32; 
flm_gen= data;      % i.e. cell2mat(flm) -> [768*1024] 
isnumeric(flm_gen)
for row = 1: rowend
   flm_data{row} =   flm_gen(row,:); 
   f_gen{row}  = ssht_inverse(flm_data{row} , L,'Reality', reality);
end


%{
for n = 1:15
flm = fitsread(filename,'binarytable',...
               'Info',info,...
               'TableColumns',[1],...
               'TableRows', [n])     

% iscell(flm)  
%whos flm
%len= length(flm)
%temp = flm{len};
%sz = size(temp)

flm_gen= cell2mat(flm);
%isnumeric(flm_gen)
L = 32;
disp('f_gen from flm_gen')
f_gen = ssht_inverse(flm_gen, L,'Reality', reality);
subplot(5, 3, n);
ssht_plot_mollweide(f_gen, L);
end 
%}
% ---------------
% Define curvelet parameters: 
% ---------------
B = 2;   
Spin = 0; 
N= L;     % Since m=l, the azimuthal band limit N = overall band limit L
J_min = 3; % minimum and maximum scale probed by wavelets 
J =s2let_jmax(L, B);  %=ceil(log L/ log B);   
 
% ---------------
% Define Plotting parameters:
% ---------------
zoomfactor =1.;  %1.6;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 16; %10  %16 % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 4;  %5 % ns;

maxfigs = nx*ny;
pltroot = '../../../figs' ; 
configstr = ['N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]; 


row = 1  % from 1 - 768
% -----------------
% Signal analysis: 
% -----------------
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_px2cur(f_gen{row},  ...
                                                                'B', B, 'L', L, ...
                                                                'J_min', J_min, ...
                                                                'Spin', Spin, ...
                                                                'Reality', reality, ...
                                                                'Upsample', upsample , ...
                                                                'SpinLowered', false, ...
                                                                'SpinLoweredFrom', 0,...
                                                                'Sampling', 'MW');


figure('Position',[20 20 1700 1400]) %100 100 1300 1000
subplot(ny, nx, 1);
% --- plot initial data --- % 
ssht_plot_mollweide(f_gen, L);
%
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])

%
subplot(ny, nx, 2);
% --- plot scaling function contributions --- %
if (upsample ~= false) %full resolution 
    ssht_plot_mollweide(f_scal, L);
else
    band_limit = min([ s2let_bandlimit(J_min-1,J_min,B,L) L ]);
    ssht_plot_mollweide(f_scal, band_limit);
end
%
title('Scaling fct')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])


% --- plot curvelet kernel contributions --- % 
ind = 2;
for j = J_min:J  
    band_limit = min([ s2let_bandlimit(j,J_min,B,L) L ]);  
    Nj = band_limit; 
    if (reality == false) % default i.e. complex signals
        enmax = 2*Nj-1; 
    else % i.e. real signals
        enmax = Nj; 
    end 
	for en = 1: enmax
		ind = ind + 1;
        if ind <= maxfigs
            subplot(ny, nx, ind);
            %
            if (upsample ~= false)
                ssht_plot_mollweide(reshape(f_cur{j-J_min+1}(en,:,:), L, 2*L-1), L);
            else
                ssht_plot_mollweide(reshape(f_cur{j-J_min+1}(en,:), band_limit, 2*band_limit-1), band_limit);
            end
            %
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Curvelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)],'FontSize', 10)
        end
    end
end
colormap(jet)
if (upsample ~= false)
    fname = [pltroot,'/s2let_demo9_', configstr, '_spin_curvelet_EarthTomo_FullRes.png']
else
    fname = [pltroot,'/s2let_demo9_', configstr, '_spin_curvelet_EarthTomo_MultiRes.png']
end
print('-r200', '-dpng', fname)

% ---------- 
% Compare reconstructed signal with the initial signals: 
% ---------- 
%{
f_rec = s2let_transform_spin_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                      'B', B, 'L', L, ...
                                                      'J_min', J_min, ...
                                                      'Spin', Spin, ...
                                                      'Reality', reality, ...
                                                      'Upsample', upsample, ...
                                                      'SpinLowered', false, ...
                                                      'SpinLoweredFrom', 0,...
                                                 	  'Sampling', 'MW');

figure('Position',[100 100 900 200]) 
subplot(2, 2, 1);
ssht_plot_mollweide(f_gen, L);
title('initial signal')
hold on
subplot(2, 2, 2);
ssht_plot_mollweide(f_rec,L);
title('reconstructed signal')
if (upsample == false)
    fname = [pltroot,'/s2let_demo9_', configstr, '_spin_curvelet_EarthTomo_FullRes_Int_Rec_signal.png']
else
    fname = [pltroot,'/s2let_demo9_', configstr, '_spin_curvelet_EarthTomo_MultiRes_Int_Rec_signal.png']
end
print('-r200', '-dpng', fname)
% Check error:
check_error = max(abs(f_gen(:)-f_rec(:)))
%}                           

% ============
% Directional wavelets
% ============
%{
[f_wav, f_scal] = s2let_transform_analysis_mw(f_gen, 'B', B, 'J_min', J_min, ...
                                              'N', N, 'Reality',reality,...
                                              'Upsample', upsample, 'Spin', 0);

% FULL-RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(ny, nx, 1);
ssht_plot_mollweide(f_gen, L);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(ny, nx, 2);
ssht_plot_mollweide(f_scal, L);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
title('Scaling fct')
ind = 2;
for j = J_min:J
	for en = 1:N
		ind = ind + 1;
        if ind <= maxfigs
            subplot(ny, nx, ind);
            ssht_plot_mollweide(f_wav{j-J_min+1,en}, L, 'Mode', 1);
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Wavelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)])
        end
	end
end
colormap(jet)
fname = [pltroot,'/s2let_demo9_', configstr, '_directional_wavelets_EarthTomo_fullres.png']
print('-r200', '-dpng', fname)
%}