function [ridgelet_wav, ridgelet_scal] = ...
   s2let_ridgelet_compute_wav(L, varargin)
% s2let_ridgelet_compute_wav - Compute ridgelets
%
% Compute ridgelets and scaling functions.
%
% Default usage :
%
%   [ridgelet_wav, ridgelet_scal] = s2let_ridgelet_compute_wav(L, <options>)
%
% where L is harmonic band-limit for the reconstruction on the sphere,
% [ridgelet_wav, ridgelet_scal] are the reconstruct ridgelet wavelets and
% scaling functions on the sphere for each j.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'B'               = { Dilation factor; B > 1 (default = 2) }
%  'Spin'            = { Spin number; Spin >= 0 (default = 0) }
%  'J_min'           = { Minimum wavelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012-2015-2014  Boris Leistedt, Martin Büttner & Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('L', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Reality', false, @islogical);
p.parse(L, varargin{:});
args = p.Results;

N = 1;
s = args.Spin; 

% Compute ring in harmonic space.
ring_lm = zeros(args.L^2,1);
for el = max([0 abs(args.Spin)]):args.L-1
        
   % Spin zero setting.
%    elon2 = el./2.0;
%    logp = gammaln(2*elon2+1) - 2*elon2 * log(2) - 2 * gammaln(elon2+1);
%    p0 = real((-1).^elon2) .* exp(logp);
   
   % Spin setting but computed via associated Legendre functions (slow).
%    P = legendre(el, 0);   
%    p0 = P(abs(args.Spin)+1);
   
   % Spin setting but computed explicitly.
   
   logp2 = gammaln(el+s+1) - el * log(2) - gammaln((el+s)./2+1) - gammaln((el-s)./2+1);
   p0 = real((-1).^((el+s)./2)) .* exp(logp2);
   
   ind = ssht_elm2ind(el, 0);
   ring_lm(ind) = 2*pi * sqrt((2*el+1)/(4*pi)) * p0;      
   ring_lm(ind) = ring_lm(ind) .* ...   
     (-1).^args.Spin .* sqrt(exp(gammaln(el-s+1) - gammaln(el+s+1)));
   
end


%% Compute ridgets via wavelet transforms
% Can only do this for scalar setting 

% % Compute ring on the sphere.
% ring = ssht_inverse(ring_lm, L, ...
%    'Reality', args.Reality, ...
%    'Method', args.Sampling, ...
%    'Spin', args.Spin);
% 
% % Compute ridgelets.
% [ridgelet_wav, ridgelet_scal] = ...
%    s2let_transform_analysis_mw(ring, 'B', args.B, 'J_min', args.J_min, ...
%    'N', N, 'Upsample', true, 'Spin', args.Spin, ...
%    'Reality', args.Reality, 'Sampling', args.Sampling);

% [ridgelet_wav, ridgelet_scal] = ...
%    s2let_transform_axisym_analysis_mw(ring, 'L', L, 'B', args.B, 'J_min', args.J_min, ...
%    'Upsample', true, ...
%    'Reality', args.Reality);


%% Compute ridgelets directly 
% Can then compute explicitly for spin functions. 
% But not that spin ridgelets are identical to scalar wavelets since can
% always view as convolution of ring with scalar wavelets.

J = s2let_jmax(args.L, args.B);
[kappa, kappa0] = s2let_transform_axisym_tiling(args.B, args.L, args.J_min);

ridgelet_scal_lm = zeros(args.L^2,1);
for el = 0:args.L-1
   ind = ssht_elm2ind(el, 0);
   scal_lm = kappa0(el+1) .* sqrt((2*el+1)/(4*pi));
   ridgelet_scal_lm(ind) = scal_lm .* sqrt(4*pi / (2*el+1)) .* ring_lm(ind);
end
ridgelet_scal = ssht_inverse(ridgelet_scal_lm, L, 'Spin', args.Spin, ...
   'Reality', args.Reality, 'Method', args.Sampling);

ridgelet_wav = cell(J-args.J_min+1,1);
for j = args.J_min:J
      
   ridgelet_wav_lm = zeros(args.L^2,1);
   for el = 0:args.L-1
      ind = ssht_elm2ind(el, 0);
      wav_lm = kappa(j+1, el+1) .* sqrt((2*el+1)/(8*pi.^2)); 
      ridgelet_wav_lm(ind) = wav_lm .* sqrt(4*pi / (2*el+1)) .* ring_lm(ind);
   end
   ridgelet_wav{j-args.J_min+1} = ssht_inverse(ridgelet_wav_lm, L, ...
      'Spin', args.Spin, 'Reality', args.Reality, 'Method', args.Sampling);
   
end






