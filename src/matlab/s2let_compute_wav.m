function psi_j = s2let_compute_wav(j, alpha, beta, gamma, L, varargin)
% s2let_compute_wav - Compute a rotated wavelet
%
% Compute the j-th wavelet, rotated by rho=(alpha, beta, gamma) in
% harmonic space and reconstruct it on the sphere.
%
% Default usage :
%
%   psi_j = s2let_compute_wav(j, alpha, beta, gamma, L, <options>)
%
% j is the order of the wavelet under consideration (depends on B)
% rho=(alpha, beta, gamma) is the rotation in SO(3) by which to rotate
% the wavelet wavelet
% L if harmonic band-limit for the reconstruction on the sphere
% psi_j is the reconstructed wavelet on the sphere, at resolution L
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'B'               = { Dilation factor; B > 1 (default = 2) }
%  'N'               = { Azimuthal band-limit; N > 0 (default = L) }
%  'Spin'            = { Spin number; Spin >= 0 (default = 0) }
%
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012-2015-2014  Boris Leistedt, Martin Büttner & Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('j', @isnumeric);
p.addRequired('alpha', @isnumeric);
p.addRequired('beta', @isnumeric);
p.addRequired('gamma', @isnumeric);
p.addRequired('L', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('N', -1, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);

p.parse(j, alpha, beta, gamma, L, varargin{:});

args = p.Results;

if args.N == -1
    args.N = L;
end

B = args.B;
N = args.N;
Spin = args.Spin;

psi = s2let_wavelet_tiling(B, L, N, Spin, j);

% Precompute Wigner small-d functions
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end

% Rotate spherical harmonic
psi_lm_rot = ssht_rotate_flm(psi(:,j+1), d, alpha, gamma);

psi_j = ssht_inverse(complex(real(psi_lm_rot), imag(psi_lm_rot)), L, 'Spin', Spin);

end
