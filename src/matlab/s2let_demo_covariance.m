% s2let_demo_covariance - Run covariance demo.
%
% Demo to compare theoretical covariance of wavelet coefficients with
% empirical data from using our transform functions. The empirical
% covariance is computed for several sets of harmonic coefficients,
% and the theoretical covariance is compared to the mean of those
% calculations in units of its standard deviation.
%
% Default usage is given by
%
%   s2let_demo_covariance
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% S2LET package to perform Wavelet transforms on the sphere.
% Copyright (C) 2012-2015-2014 Boris Leistedt and Jason McEwen
% See LICENSE.txt for license details

clear all;

% Define parameters
B = 2;
L = 16;
N = L;
Spin = 0;
J_min = 2;
J = s2let_jmax(L, B);

var_flm = 1; % Should we use the actual variance var(flmn) of each
             % realization here instead?

% Compute theoretical covariance of wavelet coefficients.
% The covariance <Wj(rho)Wj*(rho)> is given by the following expression:
%
% sigma^2 Sum(l,n) |Psij_ln|^2
%
% Where sigma^2 is the variance of the harmonic coefficients and Psij_lm
% are the harmonic coefficients of the j-th wavelet.
%
% A similar expression applies for the scaling function coefficients.

[psi_lm, phi_l] = s2let_wavelet_tiling(B, L, N, Spin, J_min);

covar_W_phi_theory = 0;
for el = 0:L-1
    covar_W_phi_theory = covar_W_phi_theory + phi_l(el+1) * phi_l(el+1)';
end

covar_W_psi_theory = zeros(J-J_min+1,1);
for j = J_min:J
    ind = 1;
    for el = 0:L-1
        for n = -el:el
            covar_W_psi_theory(j-J_min+1) = covar_W_psi_theory(j-J_min+1) + psi_lm(ind,j+1) * psi_lm(ind,j+1)';
            ind = ind + 1;
        end
    end
end

covar_W_phi_theory = covar_W_phi_theory .* var_flm
covar_W_psi_theory = covar_W_psi_theory .* var_flm

runs = 100;
if J_min == 0
    % Ergodicity fails for J_min = 0, because the scaling function will
    % only cover f00. Hence <flm flm*> will be 0 in that case and the
    % scaling coefficients will all be the same. So, if we do have
    % J_min = 0, we take the variance over all realisations instead (of
    % course, we then won't have a standard deviation to compare it to the
    % theoretical variance).
    f_scals = zeros(runs,1);
else
    covar_W_phi_data = zeros(1, runs);
end

covar_W_psi_data = zeros(J-J_min+1, runs);


for i = 1:runs
    % Generate normally distributed random flmn of complex signal
    % with mean 0 and variance 1
    flm = zeros(L^2 - Spin^2, 1);
    flm = [zeros(Spin^2,1); (randn(size(flm)) + 1i*randn(size(flm)))/sqrt(2)*sqrt(var_flm)];

    % Compute inverse then forward transform.
    f = ssht_inverse(flm, L, 'Spin', Spin);

    [f_wav, f_scal] = s2let_transform_analysis_mw(complex(real(f), imag(f)), 'B', B, 'L', L, 'N', N,    ...
                                                  'Spin', Spin, 'J_min', J_min, ...
                                                  'Upsample', true);

    if J_min == 0
        f_scals(i) = f_scal(1);
    else
        covar_W_phi_data(i) = var(f_scal(:));
    end

    for j = J_min:J
        f_wav_j = f_wav{j-J_min+1};
        covar_W_psi_data(j-J_min+1,i) = var(f_wav_j(:));
    end
end

if J_min == 0
    covar_W_phi_data = var(f_scals(:))
    W_phi_error_absolute = abs(covar_W_phi_data - covar_W_phi_theory)
else
    mean_covar_W_phi_data = mean(covar_W_phi_data)
    std_covar_W_phi_data = std(covar_W_phi_data)
    W_phi_error_in_std = abs(mean_covar_W_phi_data - covar_W_phi_theory)/std_covar_W_phi_data
end

mean_covar_W_psi_data = mean(covar_W_psi_data,2)
std_covar_W_psi_data = std(covar_W_psi_data,0,2)
W_psi_error_in_std = abs(mean_covar_W_psi_data - covar_W_psi_theory)./std_covar_W_psi_data
