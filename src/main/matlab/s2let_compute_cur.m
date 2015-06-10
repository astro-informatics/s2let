function cur_lm_j = s2let_compute_cur(j, alpha, beta, gamma, L, varargin)
% s2let_compute_cur - Compute a rotated curvelet
%
% Compute the j-th curvelet, rotated by rho=(alpha, beta, gamma) in
% harmonic space and reconstruct it on the sphere.
%
% Default usage :
%
%   cur_lm_j = s2let_compute_cur(j, alpha, beta, gamma, L, <options>)
%
% j is the order of the curvelet under consideration (depends on B)
% rho=(alpha, beta, gamma) is the rotation in SO(3) by which to rotate
% the curvelet
% L if harmonic band-limit for the reconstruction on the sphere
% psi_j is the reconstructed curvelet on the sphere, at resolution L
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'B'               = { Dilation factor; B > 1 (default = 2) }
%  'N'               = { Azimuthal band-limit; N > 0 (default = L) }
%  'Spin'            = { Spin number; Spin >= 0 (default = 0) }
%
% S2LET package to perform curvelet transform on the Sphere.
% Copyright (C) 2012-2014  Boris Leistedt, Martin Büttner & Jason McEwen
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

cur_lm = s2let_curvelet_tiling(B, L, j, 'Spin', Spin);

% Precompute Wigner small-d functions
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end

% Rotate spherical harmonic
cur_lm_rot = ssht_rotate_flm(cur_lm{j-J_min+1}(:), d, alpha, gamma);

cur_lm_j = ssht_inverse(complex(real(cur_lm_rot), imag(cur_lm_rot)), L, 'Spin', Spin);

end
