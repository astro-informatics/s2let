% s2let_curvelet_matlab_prototype
% Run curvelet analysis and synthesis 
% of randomly generated signals f 

clear all;
close all;

% Main parameters
Spin = 1;
L = 16;
N = L;  % For curvelets: N=L, %N=1 for axisymmetric curvelets
B = 2;
J_min = 0;
J =s2let_jmax(L, B);  %=ceil(log L/ log B);  

%%disp('Generates random band-limited function')
%flm_gen = zeros(L^2,1);
%flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
%flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);
%disp('Construct the corresponding signal on the sphere')
%%** using C: 
%f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');


%open file
fid= fopen('/Users/jenniferyhchan/WaveletsCode_PhD/s2let/1_flm_randgen_mw_test.dat');
rawData=fscanf(fid, '%f, %f',[2 256]);
%fclose(fid);
complexData=complex(rawData(1,:),rawData(2,:));
% Non-conjugate transpose
complexData= complexData.' ; 

% Read from file "1_cur_flm_randgen_mw_test.dat" containing randomly-generated flm
flm_gen = zeros(L^2,1);
flm_gen =  complexData ; 

f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');


%disp('Perform spin directional harmonic-to-curvelet (lm2cur) transform with custom parameters')
%[f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
%flm_rec = s2let_transform_synthesis_lm2cur(f_cur, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
%default = max(abs(flm-flm_rec))

%***** BUILD PROPTOTYPE of s2let_transform_analysis_lm2cur
%      (i.e. s2let_analysis_lm2cur in s2let_analysis.c) ****************
%***** 1st) s2let_tiling_curvelet(cur_lm, scal_l, parameters);
%***** step 1a) call curvelet- and scaling-function- generating functions 
disp('Tile curvelets in harmonic space (cur_lm, scal_l)')
%[kappa kappa0] =  s2let_transform_axisym_tiling(B, L, J_min);
%[psi_lm phi_l] = s2let_wavelet_tiling(B, L, N, Spin, J_min);
[psi_lm phi_l] = s2let_curvelet_tiling(B, L, N, Spin, J_min);

%% Plot
%% Rotate function on the sphere such that x-axis -> z-axis  (beta =pi/2) 
%% (i.e. curvelets lie on the sphere top)
alpha =  pi 
beta = pi/2 
gamma = 0
%% Precompute Wigner small-d functions
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end
zoomfactor = 1.4;
plot_caxis_scale = 2
type = 'colour';
lighting = true;
nx = 3;
ny = 4;
maxfigs = nx*ny;
%% 

%% Plot curvelets on the sphere
figure('Position',[100 100 1200 600])
%cur_lm = zeros(L^2,1);
ind_pm = 0;
ind_nm = 0;
ind=0;    %for plotting
for j = J_min:J, 
 cur_lm{j-J_min+1} = psi_lm(:,j+1);  
 %cur_lm{j-J_min+1} = zeros(L^2,1);
%for el = 0:L-1, 
 %  m = el;
  %%for positive m 
 %  ind_pm = ssht_elm2ind(el, m);
  %% Curvelet coefficients: 
 %  cur_lm{j-J_min+1}(ind_pm) = kappa(j+1,el+1);  
  %%for negative m 
 %  ind_nm = ssht_elm2ind(el, -m);
 %  cur_lm{j-J_min+1}(ind_nm) = (-1)^m * conj(cur_lm{j-J_min+1}(ind_pm));
%end

%% Rotate the curvelets coefficients 
   flm_cur_rot = ssht_rotate_flm(cur_lm{j-J_min+1}, d, alpha, gamma);
   if Spin == 0
%% Compute the function (rotated): 
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
%***** step 1b) compute the scaling coefficients (no j-dependence except J_min)
scal_l = zeros(L^2,1);
for l = 0:L-1, 
%scal_l(l^2+l+1,1) = kappa0(l+1);
scal_l(l^2+l+1,1) = phi_l(l+1);
end
%%%  f_scal = ssht_inverse(scal_l, L, 'Reality', true);

%**** 2nd) s2let_analysis_lm2lmn(f_cur_lmn, f_scal_lm, flm, cur_lm, scal_l, parameters);
% Generate flmn of complex signal 
disp(' ')
disp('Signal Analysis: ')
disp('Curvelet analysis of complex signals in Wigner space (i.e. flm to flmn)')
%zeros((2*N-1)*(L*L-N*(N-1)/3), 1);
%(see  so3_core_inverse_via_ssht.m for the size)
%offset = 0;
ind_ln=0; 
ind = 0;
ind_lmn = 0;
for j = J_min:J,
 flmn{j-J_min+1} = zeros((2*N-1)*L*L,1);
 for n = -N+1:N-1,
    for el = abs(n):L-1,
        ind_ln = ssht_elm2ind(el, n);   %for l and n 
        psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-J_min+1}(ind_ln));  % 
        for m = -el:el,
        ind = ssht_elm2ind(el, m); 
        ind_lmn = so3_elmn2ind(el,m,n,L,N);
        flmn{j-J_min+1}(ind_lmn) =  flm_gen(ind) * psi; 
%        flmn(offset+ind) =  flm_gen(ind) * psi;        
        end
    end
 end
% offset =  offset+ (2*N-1)*L*L ;    % i.e. for (padded, complex)
%! Remark about offset in C-code: offset += so3_sampling_flmn_size (L, N)
%! C code: so3_sampling_flmn_size (L, N)
%! Padded :  if Real: N*L*L; if Complex: (2*N-1)*L*L
%! Compact:  if Real: N*(6*L*L-(N-1)*(2*N-1))/6; if Complex:(2*N-1)*(3*L*L-N*(N-1))/3;
%! Both of these are based on the fact that the sum
%! over n*n from 1 to N-1 is (N-1)*N*(2*N-1)/6.%so3_sampling_flmn_size(&so3_parameters);
end
%% Scaling function
disp('Compute scaling function f_scal_lm=flm_gen(lm_ind) * phi ')
%f_scal_lm = zeros(L^2,1)
lm_ind=0;
for j = J_min:J, 
 for el = 0:L-1,   
  %phi = sqrt(4.0*pi/(2.*el+1)) * kappa0(el+1);
  phi = sqrt(4.0*pi/(2.*el+1)) * phi_l(el+1);  %
  for m = -el:el, 
  lm_ind=ssht_elm2ind(el, m);
  f_scal_lm(lm_ind) = flm_gen(lm_ind) * phi;
  end
 end
end
    
% Compute inverse then forward transform.
disp('so3_forward (i.e. f to flmn) and so3_inverse (i.e. flmn to f): ')
for j = J_min:J, 
f_cur = so3_inverse(flmn{j-J_min+1}, L, N);  %i.e. so3_core_inverse_via_ssht -> f_cur[nalpha*nbeta*ngamma*sizeof(f)]
% f_cur+offset = so3_inverse(flmn+offset_lmn, L, N);
% where in C code: offset_lmn += so3_sampling_flmn_size
%                  offset += so3_sampling_f_size
flmn_syn{j-J_min+1} = so3_forward(f_cur, L, N);%i.e. so3_core_forward_via_ssht -> input funcion on the sphere and output flmn harmonic coefficient 
% Compute maximum error in harmonic space.
disp('Check the error of so3_forward (i.e. f to flmn) and so3_inverse (i.e. flmn to f): ')
maxerr = max(abs(flmn_syn{j-J_min+1} - flmn{j-J_min+1}))
end

% void s2let_synthesis_cur_lmn2lm(
disp(' ')
disp('Signal synthesis: ')
% Compute flm_syn
disp('Compute flm_cur_syn from flmn_syn')
%offset = 0;
ind_ln =0; 
ind=0; 
ind_lmn=0; 
flm_cur_syn=zeros(L^2,1);
for j = J_min:J, 
 for n = -N+1:N-1,
    for el = abs(n):L-1,
        ind_ln = ssht_elm2ind(el, n);   %for l and n 
        psi = (cur_lm{j-J_min+1}(ind_ln));  % no ((2*el+1)/(8.*pi*pi)) in C
        for m = -el:el,
        ind = ssht_elm2ind(el, m);   %for positive m
        ind_lmn = so3_elmn2ind(el,m,n,L,N);
%        flm_test_cur(ind) =flm_test_cur(ind)+ flmn(offset+ind_lmn)* psi;   % sum over m (c.f. n)
        flm_cur_syn(ind) =flm_cur_syn(ind)+ flmn{j-J_min+1}(ind_lmn)* psi;   % sum over m (c.f. n)
        end 
    end
 end
% offset =  offset+ (2*N-1)*L*L ;    % i.e. for (padded, complex)
end
disp('Compute flm_scal_syn ')
%% For scaling function - sum over m: 
flm_scal_syn=zeros(L^2,1);
lm_ind=0;  
%for j = J_min:J, 
 for el = 0:L-1,
 % phi = sqrt(4.*pi/(2.*el+1)) * kappa0(el+1);
   phi = sqrt(4.*pi/(2.*el+1)) * phi_l(el+1);  %
  for m = -el:el, 
  lm_ind=ssht_elm2ind(el, m);
  flm_scal_syn(lm_ind) =  flm_scal_syn(lm_ind)+ flm_gen(lm_ind) * phi;
  end
 end
%end
% f_scal_l=sum(flm_test_scal)


disp('Summing flm_cur_syn+flm_scal_syn ')
disp('then ')
disp('Compute the re-constructed function via ssht_inverse ')

 flm_test=flm_scal_syn+flm_cur_syn;
 f_test_syn = ssht_inverse(flm_test, L, 'Method', 'MW'); 


disp('Done- check results ')
disp('Check the difference between flm_gen and flm_test: ')
maxerr = max(abs(flm_gen - flm_test))
%disp('Check the difference between f_gen and f_test_syn: ')
%maxerr = max(abs(f_gen - f_test_syn))
