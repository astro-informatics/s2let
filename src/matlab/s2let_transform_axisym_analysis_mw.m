function [f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f, varargin)

% s2let_transform_axisym_analysis_mw
% Compute axisymmetric wavelet transform, output in pixel space.
%
% Default usage :
%
%   [f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f, <options>)
%
% f is the input field -- MW sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'J_min'           = { Minimum wavelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true       [full resolution wavelets] }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
Lguessed = min([sz(1) sz(2)]);

p = inputParser;
p.addRequired('f', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Reality', false, @islogical);
p.parse(f, varargin{:});
args = p.Results;

f_vec = s2let_mw_arr2vec(f);

[f_wav_vec, f_scal_vec] = s2let_transform_axisym_analysis_mw_mex(f_vec, args.B, args.L, args.J_min, args.Reality, args.Upsample);

f_scal = s2let_mw_vec2arr(f_scal_vec);

J = s2let_jmax(args.L, args.B);
f_wav = cell(J+1-args.J_min, 1);
offset = 0;
for j = args.J_min:J
  if args.Upsample
    band_limit = args.L;
  else
    band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
  end
  temp = zeros(band_limit, 2*band_limit-1);
  for t = 0:band_limit-1
      for p = 0:2*band_limit-2
          ind = offset + t * ( 2 * band_limit - 1) + p + 1;
          temp(t+1,p+1) = f_wav_vec(1,ind);
      end
  end
  f_wav{j+1-args.J_min} = temp;
  offset = offset + band_limit * (2*band_limit-1);
end
