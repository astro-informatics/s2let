% s2let_curvelet_matlab_prototype
% Run curvelet analysis and synthesis 
% of randomly generated signals f 

clear all;
close all;

% Main parameters
Spin = 0;
L = 16;
N = L; % For curvelets: N=L  %N=1 for axisymmetric curvelets
B = 2;
J_min = 1;  %N.B. J_min = 2 : NOT EXACT
J =s2let_jmax(L, B);  %=ceil(log L/ log B);  

%{
disp('Generates random band-limited function')
flm_gen = zeros(L^2,1);
flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);

disp('Construct the corresponding signal on the sphere')
f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');
%}

flm_gen = zeros(L^2,1);
disp('read from file flm_gen'); 
fid= fopen('/Users/jenniferyhchan/WaveletsCode_PhD/s2let_curvelets_MATLAB/1_cur_flm_randgen_mw_test.dat');
rawData=fscanf(fid, '%f, %f',[2 256]);
fclose(fid);
complexData=complex(rawData(1,:),rawData(2,:));
flm_gen= complexData.' ;  % Non-conjugate transpose
disp('Construct the corresponding signal via ssht_inverse'); 
f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');
flm_gen= ssht_forward(f_gen, L, 'Method', 'MW');

%
% The following c.f.the C analysis function:
%  s2let_analysis_px2wav( complex double *f_scal, const complex double *f,
%                          const s2let_parameters_t *parameters) 
%
% * Wavelet analysis from pixel space to wavelet space for complex signals.
% *
% * \param[out]  f_wav Array of wavelet maps
% * \param[out]  f_scal Scaling function map
% * \param[in]   f Signal on the sphere
% * \param[in]   parameters A fully populated parameters object.

% ---------------
% Tile curvelets:
% ---------------
%***** step 1a) call curvelet- and scaling-function- generating functions 
disp('Tile curvelets in harmonic space (cur_lm, scal_l)')
[psi_lm phi_l] = s2let_curvelet_tiling(B, L, N, Spin, J_min);
for j = J_min:J, 
 cur_lm{j-J_min+1} = psi_lm(:,j+1);  
end
%***** step 1b) compute the scaling coefficients (no j-dependence except on J_min)
scal_l = zeros(L^2,1);
for l = 0:L-1,
scal_l(l^2+l+1,1) = phi_l(l+1);
end

% -----------------
% Signal analysis: (harmonic to pixel space) 
% -----------------
% Call matlab function analysis_lm2cur
[f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm_gen,  ...
                                                  'B', B, 'L', L, ...
                                                  'Spin', Spin, 'J_min', J_min, ...
                                                  'N', N, 'Sampling', 'MW',...
                                                  'Upsample', false);
                                              
% -----------------
% Signal synthesis: (pixel to harmonic space) 
% -----------------
% Call s2let_transform_synthesis_lmn2lm
%{
[flm_cur_syn, flm_scal_syn] = s2let_transform_synthesis_lm2cur(f_cur, f_scal, flm_gen,...
                                           'B', B, 'L', L, 'J_min', J_min, ...
                                           'N', N,'Sampling', 'MW', ...
                                           'Upsample', false, 'Spin', 0);
disp('Sum: flm_cur_syn+flm_scal_syn ');
flm_rec = flm_scal_syn+flm_cur_syn;
                                       
disp('Compute the re-constructed function via ssht_inverse ');
f_rec = ssht_inverse(flm_rec, L, 'Method', 'MW'); 

disp('Both analysis and synthesis are done! ');
disp('');


% ---- Check whether tranform is exact --- % 
disp('- Test exact transform: check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))
%}

% ------------------- 
% FULL RESOLUTION PLOT
% ------------------- 

zoomfactor = 1.6;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 4; % ns - 1 + rem(2+J-J_min+1 , ns) ;
nx = 3; % ns;

maxfigs = nx*ny;
pltroot = '../../../figs'
configstr = ['N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]

figure('Position',[100 100 1300 1000])
subplot(ny, nx, 1);
ssht_plot_mollweide(f_gen, L, 'Mode', 1);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
subplot(ny, nx, 2);
ssht_plot_mollweide(f_scal, L, 'Mode', 1);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
title('Scaling fct')
ind = 2
for j = J_min:J
	for en = 1:2*N-1
		ind = ind + 1
        if ind <= maxfigs
            subplot(ny, nx, ind);
            ssht_plot_mollweide(f_cur{j-J_min+1,en}, L, 'Mode', 1);
            campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
            v = caxis;
            temp = max(abs(v));
            caxis([-temp temp])
            title(['Wavelet scale j=',int2str(j)-J_min+1,', n=',int2str(en)])
        end
	end
end

colormap(jet)
fname = [pltroot,'/s2let_demo4_', configstr, '_earth_multires.png']
print('-r200', '-dpng', fname)


