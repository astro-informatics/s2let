% s2let_proptotype_fulltest
% Run all exactness tests for the MW sampling,
% all wavelet transforms must reconstruct the input maps
% at floating-point precision. Various parameters are tested.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 16;
N = L;
B = 2;
Spin = 0;
J_min = 2;
J =s2let_jmax(L, B);  %=ceil(log L/ log B);  

disp('Generates random band-limited function')
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
%** using C: 
%** f = ssht_inverse(flm, L, 'Method', 'MW');
% using matlab prototype "inverse_direct" , where
% function [ f ] = inverse_direct( L, flm )
% INVERSE_DIRECT Inverse spherical harmonic transform using only matlab and
% a simple loop.
%.
%   This computes the inverse spherical harmonic transform, simply
%   synthesising the function by summing over all spherical harmonics.
%
%   L   ... the band limit (maximum l is L-1)
%   flm ... the coefficients of the spherical harmonics
%           to avoid wasting memory, these should be supplied in an
%           unrolled array of the following format:
%           [(0,0) (1,-1) (1,0) (1,1) (2,-2) (2,-1) ... ]
%           of size L^2, where the first number corresponds to l
%           and the second to m.

if length(flm) ~= L^2
    error('Parameter flm has to contain L^2 coefficients.')
end

[thetas, phis] = ssht_sampling(L);
f_gen = zeros(length(thetas), length(phis));
for j=1:length(thetas),
    theta = thetas(j);
    for k=1:length(phis),
        phi = phis(k);
        sum = 0;
        lmindex = 1;
        for l = 0:L-1,
            ppos = legendre(l,cos(theta));
            for m = -l:l,
                if m >= 0
                    p = ppos(m+1);
                else
                    p = (-1)^m*factorial(l+m)/factorial(l-m)*ppos(-m+1);
                end
                sum = sum + flm(lmindex)*...
                    sqrt(...
                      (2*l+1)*factorial(l-m)/...
                      ((4*pi)*factorial(l+m))...
                    )*...
                    p*exp(1i*m*phi);
                lmindex = lmindex + 1;
            end
        end
        f_gen(j,k) = sum;
    end
end

%disp('Perform spin directional harmonic-to-wavelet (lm2wav) transform with custom parameters')
%[f_wav, f_scal] = s2let_transform_analysis_lm2wav(flm,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
%flm_rec = s2let_transform_synthesis_lm2wav(f_wav, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
%default = max(abs(flm-flm_rec))

%***** BUILD PROPTOTYPE of s2let_transform_analysis_lm2wav
%      (i.e. s2let_analysis_lm2wav in s2let_analysis.c) ****************
%***** 1st) s2let_tiling_wavelet(wav_lm, scal_l, parameters);
%***** step 1a) call wavelet- and scaling-function- generating functions 
disp('Tile curvelets in harmonic space (wav_lm, scal_l)')
[kappa kappa0] =  s2let_transform_axisym_tiling(B, L, J_min);
% step 1b) compute the curvelet coefficients in harmonic space: wav_lm
wav_lm = zeros(L^2,1);
ind_pm = 0;
ind_nm = 0;
for j = J_min:J
for el = 0:L-1
m = el;
ind_pm = ssht_elm2ind(el, m);  %for positive m
disp(ind_pm)
% Curvelet coefficients:
wav_lm(ind_pm) = kappa(j+1,el+1);
ind_nm = ssht_elm2ind(el, -m);  %for negative m
wav_lm(ind_nm) = (-1)^m * conj(wav_lm(ind_pm));
end
end
%  Compute the wavelet function :
%f_wav= ssht_inverse(wav_lm, L, 'Reality', false);


%% Rotate function on the sphere such that x-axis -> z-axis  (beta =pi/2) 
% (i.e. curvelets lie on the sphere top)
alpha = pi 
beta = pi/2 
gamma = 0
% Precompute Wigner small-d functions
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end

%% Plotting: 
% Plot curvelets on the sphere
figure('Position',[100 100 1200 600])
zoomfactor = 1.4;
plot_caxis_scale = 2
type = 'colour';
lighting = true;
%ns = ceil(sqrt(2+J-J_min+1)) ;
%nx = ns - 1 + rem(2+J-J_min + 1, ns) ;
%ny = ns;
nx = 3
ny = 4

maxfigs = nx*ny;
% Rotate the wavelets coefficients 
flm_wav_rot = ssht_rotate_flm(wav_lm, d, alpha, gamma);
ind=0;
for j = J_min:J
   if Spin == 0
% Compute the function (rotated): 
       f_wav_rot = ssht_inverse(flm_wav_rot, L, 'Reality', true);
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
% Plot the non-rotated function:            
%           ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);

% Plot the rotated function on the sphere
           ssht_plot_sphere(f_wav_rot, L, 'Type', type, 'Lighting', lighting);
           title(h, ['Curvelet j = ',int2str(j-J_min+1)])
           locate = get(h,'title');
           pos = get(locate,'position');
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           v = caxis;
           temp = max(abs(v));
           caxis([-temp temp]*plot_caxis_scale)
           zoom(zoomfactor)
% view along the x-axis 
%           view([-90,0])
       end
   end
end

%colormap(jet)
%%fname = [pltroot,'/s2let_demo7_', configstr, '_cur_jet.png']
%fname = ['s2let_demo7_', configstr, '_cur_jet.png']
%print('-r200', '-dpng', fname)




%***** step 1b) compute the scaling coefficients
scal_l = zeros(L^2,1);
for l = 0:L-1
scal_l(l^2+l+1,1) = kappa0(l+1);
end

%**** 2nd) s2let_analysis_lm2lmn(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, parameters);
%    where in the C language: 
%    int offset = 0;
%    for (j = J_min; j <= J; ++j)
%    {
%        if (!parameters->upsample)
%        {
%            bandlimit = MIN(s2let_bandlimit(j, parameters), L);
%            so3_parameters.L = bandlimit;
%            Nj = MIN(N,bandlimit);
%            // ensure N and Nj are both even or both odd
%            Nj += (Nj+N)%2;
%            so3_parameters.N = Nj;
%        }
%        for (n = -Nj+1; n < Nj; n+=1)  //original n+=2
%{
%    for (el = MAX(ABS(spin), ABS(n)); el < bandlimit; ++el)
%    {
%        ssht_sampling_elm2ind(&lm_ind, el, n);
%        psi = 8*PI*PI/(2*el+1) * conj(wav_lm[j*L*L + lm_ind]);
%        for (m = -el; m <= el; ++m)
%        {
%            ssht_sampling_elm2ind(&lm_ind, el, m);
%            so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
%            f_wav_lmn[offset + lmn_ind] = flm[lm_ind] * psi;
%        }
%    }
%}
% offset += so3_sampling_flmn_size(&so3_parameters);
%}
%***** and scaling function :
%for (el = ABS(spin); el < bandlimit; ++el)
%{
%    phi = sqrt(4.0*PI/(2*el+1)) * scal_l[el];
%    for (m = -el; m <= el; ++m)    %************ m=0 treatment?!?!?!
%    {
%        ssht_sampling_elm2ind(&lm_ind, el, m);
%        f_scal_lm[lm_ind] = flm[lm_ind] * phi;
%    }
%}
%
% Now in matlab: 
% Generate flmn of complex signal 
disp('Curvelet analysis of complex signals in Wigner space (i.e. flm to flmn)')
flmn = zeros((2*N-1)*(L*L-N*(N-1)/3), 1); %(see  so3_core_inverse_via_ssht.m for the size)
for j = J_min:J
 for n = -N+1:N-1,
    for el = abs(n):L-1,
        m = el;   % for m = -el:el,
        %ind = so3_elmn2ind(el,m,n,L,N);
        ind_pm = ssht_elm2ind(el, m);   %for positive m
        psi_pm = 8.*pi*pi/(2*el+1) * conj(wav_lm(ind_pm));
        flmn(ind_pm) =  flm(ind_pm) * psi_pm;
        ind_nm = ssht_elm2ind(el, -m);  %for negative m
        psi_nm = (8*pi*pi/(2*el+1)) * conj(wav_lm(ind_nm));
        flmn(ind_nm) =  flm(ind_nm) * psi_nm;
    end
 end
end



%%% For scaling function: 
% 3rd) ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, bandlimit, 0, dl_method, verbosity);
%  Compute the scaling function :
f_scal = ssht_inverse(scal_l, L, 'Reality', false);

%% Plot scaling function which has no directional (i.e. m) dependence, 
%% the same as the axisymmetric wavelets and directional wavelets 
figure('Position',[100 100 300 300])
h=subplot(1, 1, 1);
ssht_plot_sphere(f_scal, L, 'Type', type, 'Lighting',  lighting);
title(h,'Scaling fct')
locate = get(h,'title');
pos = get(locate,'position');
pos(1,2) = pos(1,2)+0.7;
pos(1,1) = pos(1,1)-0.7;
set(locate,'pos',pos);
zoom(1.2)
v = caxis;
temp = max(abs(v));
caxis([-temp temp]*plot_caxis_scale)

%colormap(jet)
%%fname = [pltroot,'/s2let_demo7_', configstr, '_scal_jet.png']
%fname = ['s2let_demo7_', configstr, '_scal_jet.png']
%print('-r200', '-dpng', fname)


%% 4th) so3_core_inverse_via_ssht(f_wav + offset,f_wav_lmn + offset_lmn,&so3_parameters);
% where offset_lmn += so3_sampling_flmn_size(&so3_parameters);
%       offset += so3_sampling_f_size(&so3_parameters);
% Compute inverse then forward transform.
% f = so3_inverse(flmn, L, N);
% In prototype file: "so3_via_ssht_inverse.m": 
% function [ f ] = so3_via_ssht_inverse( L, N, flmn )
% SO3_VIA_SSHT_INVERSE Inverse fast Wigner transform exploiting the fast
% spin spherical harmonic transform.
%
%   This computes the inverse Wigner transform of a function defined on the
%   rotation group SO(3). This implementation makes use of MATLAB's ifft()
%   functions.
%   
%   L    ... band limit in beta   (|b| < L)
%            this implicitly imposes a band limit M = L on alpha
%   N    ... band limit in gamma  (|c| < N)
%   flmn ... the coefficients of the Wigner D-function conjugates.
%            To save memory and allow for easy exploitation of the fast
%            SSHT, the coefficients are expected in an unrolled form as
%            follows:
%            The flmn array consists of 2*N-1 sections, each corresponding
%            to a value of n from -N+1 to N-1.
%            Each of those sections contains L-|n| subsections, each
%            corresponding to a value of l from |n| to L-1.
%            Each of those subsections lists 2*l+1 values of m ranging 
%            from -l+1 to l-1.
%            Note that each n-section contains L²-n² values.
%            This sums to a total of (2N-1)(L²-N(N-1)/3) values.
disp('Compute inverse Wigner transform (i.e. flmn to f)')
if L < N
    error('Band limit N must not be greater than L.')
end

if sum(size(flmn) ~= [(2*N-1)*(3*L^2-N*(N-1))/3, 1]) > 0
    error('Parameter flmn has to be a column vector containing (2N-1)(L²-N(N-1)/3) coefficients (L²-n² for each n).')
end

% step 1
% SSHT of flmn over l and m to obtain fn(a,b)
% Create a mask of scaling values that can be applied to an unrolled flm
% array.
scale = zeros(L^2,1);
offset = 1; % marks the beginning of the current l-subsection
for l = 0:L-1
    scale(offset:offset+2*l) = sqrt((2*l+1)/(16*pi^3));
    offset = offset + 2*l+1;
end

fn_inverse_so3 = zeros(2*N-1, L, 2*L-1);

offset = 1; % marks the beginning of the current n-section
for n = -N+1:N-1
    % Slice the current section from flmn and pad with zeros
    flm_inverse_so3 = [zeros(n^2,1); flmn(offset:offset+L^2-n^2-1)];
    % Scale flm
    flm_inverse_so3 = flm_inverse_so3.*scale;
    
    fn_inverse_so3(N+n,:,:) = (-1)^n*ssht_inverse(flm_inverse_so3, L, 'Spin', -n);
    offset = offset + L^2 - n^2;
end

% layout of fn_inverse_so3:
% 1st index is n, from -N+1 to N-1
% 2nd index is b, from 0 to L-1
% 3rd index is a, from 0 to 2*L-2 (or 0 to L-1, and then -L+1 to -1)

% step 2
% IFFT of fn(a,b) over n to obtain f(a,b,g) (result)
% First, we shift the 0-frequency component (n = 0) to the 
% beginning of fn (as MATLAB expects it).
% Finally, we rearrange the dimensions of f, to obtain the desired layout.
f = ifft(ifftshift(fn_inverse_so3,1))*(2*N-1);
f = permute(f,[3 2 1]);
% layout of f:
% 1st index is a, from 0 to 2*L-2 (or 0 to L-1, and then -L+1 to -1)
% 2nd index is b, from 0 to L-1
% 3rd index is g, from 0 to 2*N-2 (or 0 to N-1, and then -N+1 to -1)
% ---------------------------

% NOW consider flmn_syn = so3_forward(f, L, N);
% In prototype file: " so3_via_ssht_forward.m "
% function [ flmn ] = so3_via_ssht_forward( L, N, f )
% SO3_VIA_SSHT_FORWARD Forward fast Wigner transform exploiting the fast
% spin spherical harmonic transform.
%
%   This computes the forward Wigner transform of a function defined on the
%   rotation group SO(3). This implementation makes use of MATLAB's fft()
%   functions.
%   
%   L ... band limit in beta   (|b| < L)
%         this implicitly imposes a band limit M = L on alpha
%   N ... band limit in gamma  (|c| < N)
%   f ... sampled function values. The dimensions of f correspond to alpha,
%         beta, gamma and the values should correspond to sampling points
%         as given by so3_sampling().
%
%   To save memory and allow for easy exploitation of the fast SSHT, the 
%   returned coefficients flmn are expected in an unrolled form as follows:
%   The flmn array consists of 2*N-1 sections, each corresponding
%   to a value of n from -N+1 to N-1.
%   Each of those sections contains L-|n| subsections, each
%   corresponding to a value of l from |n| to L-1.
%   Each of those subsections lists 2*l+1 values of m ranging 
%   from -l+1 to l-1.
%   Note that each n-section contains L²-n² values.
%   This sums to a total of (2N-1)(L²-N(N-1)/3) values.
disp('Compute forward Wigner transform (i.e. f to flmn)')
if L < N
    error('Parameter N must not be greater than L.')
end
[na, nb, ng] = size(f);
if na ~= 2*L-1 || nb ~= L || ng ~= 2*N-1
    error('Parameter f has to contain (2*L-1) * L * (2*N-1) coefficients.')
end
% layout of f:
% 1st index is a, from 0 to 2*M-2 (or 0 to M-1, and then -M+1 to -1)
% 2nd index is b, from 0 to L-1
% 3rd index is g, from 0 to 2*N-2 (or 0 to N-1, and then -N+1 to -1)

% step 1
% FFT of f(a,b,g) over g to obtain fn(a,b)
fn = fftshift(fft(f,[],3),3)*2*pi/(2*N-1);
% layout of fn:
% 1st index is a, from 0 to 2*L-2 (or 0 to L-1, and then -L+1 to -1)
% 2nd index is b, from 0 to L-1
% 3rd index is n, from -N+1 to N-1

% step 2
% SSHT of fn(a,b) over a and b to obtain flmn (result)
flmn_syn = zeros((2*N-1)*(3*L^2-N*(N-1))/3,1);
% Create a mask of scaling values that can be applied to an unrolled flm
% array.
scale = zeros(L^2,1);
offset = 1; % marks the beginning of the current l-subsection
for l = 0:L-1
    scale(offset:offset+2*l) = sqrt(4*pi/(2*l+1));
    offset = offset + 2*l+1;
end

offset = 1; % marks the beginning of the current n-section
for n = -N+1:N-1
    flm_forward_so3 = (-1)^n*ssht_forward(fn(:,:,N+n).', L, 'Spin', -n);
    % Scale flm
    flm_forward_so3 = flm_forward_so3.*scale;
    
    flmn_syn(offset:offset+L^2-n^2-1) = flm_forward_so3(n^2+1:end);
    
    offset = offset + L^2 - n^2;
end





% Compute maximum error in Wigner coefficents.
disp('Check error of inverse and forward Wigner transform (i.e. compare flmn)')
maxerr = max(abs(flmn_syn - flmn))

% Compute maximum error in harmonic coefficients.
disp('Check error in harmonic space (i.e. compare flm)')
maxerr = max(abs(flm_forward_so3 - flm_inverse_so3))

%% Compute maximum error in pixel space.
%% Compute sampling grids.
%%[alphas, betas, gammas, n, nalpha, nbeta, ngamma] = so3_sampling(L, N, 'Grid', true);
%ngamma = 2*N-1
%for i = 1:ngamma,
%    f_g = squeeze(f(i,:,:));
 %%  subplot(1,ngamma,i)
 %%   ssht_plot_sphere(abs(f_g), L);
 %%   title(sprintf('g = %d; |f|', i-1))
%end
%disp('Check error of the reconstruction of signals in pixel space (i.e. compare flmn)')
%maxerr = max(abs(f_g - f_gen))

