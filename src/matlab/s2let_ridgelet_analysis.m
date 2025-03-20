
function [f_ridge_wav, f_ridge_scal] = s2let_ridgelet_analysis(f, varargin)

% s2let_ridgelet_analysis
% Compute ridgelet transform on the sphere.
%
% Default usage :
%
%   [f_ridge_wav, f_ridge_scal] = s2let_ridgelet_analysis(f, varargin)
%
% where f is the real space function on the sphere to analyse, 
% f_wav contains the output wavelet ridgelet contributions, and
% f_scal contains the output scaling ridgelet contributions,
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'N'               = { Azimuthal/directional band-limit; N > 1 (default=L) }
%  'Spin'            = { Spin; (default=0) }
%  'J_min'           = { Minimum wavelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true       [full resolution wavelets] }
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'OriginalSpin' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the wavelets
%                       should be lowered from (default = 0)]
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
if sz(1) == 2*sz(2)-1 || sz(2) == 2*sz(1)-1
    Lguessed = min([sz(1) sz(2)]);
else
    Lguessed = min([sz(1) sz(2)])-1;
end

p = inputParser;
p.addRequired('f', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('OriginalSpin', 0, @isnumeric);
p.parse(f, varargin{:});
args = p.Results;

N = 1;

f_lm = ssht_forward(f, args.L, ...
   'Method', args.Sampling, ...
   'Reality', args.Reality, 'Spin', args.Spin);

f_radon_lm = s2let_radon_transform(f_lm, ...
   'Reality', args.Reality, 'Spin', args.Spin);

f_radon = ssht_inverse(f_radon_lm, args.L, ...
   'Method', args.Sampling, ...
   'Reality', args.Reality, 'Spin', 0);

[f_ridge_wav, f_ridge_scal] = ...
   s2let_transform_analysis_mw(f_radon, 'L', args.L, 'B', args.B, ...
   'J_min', args.J_min, ...
   'N', N, 'Upsample', args.Upsample, 'Spin', 0, ...
   'Reality', args.Reality, 'Sampling', args.Sampling, ...
   'OriginalSpin', args.OriginalSpin);
