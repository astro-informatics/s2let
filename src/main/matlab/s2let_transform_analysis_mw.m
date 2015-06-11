function [f_wav, f_scal] = s2let_transform_analysis_mw(f, varargin)

% s2let_transform_analysis_mw
% Compute spin directional wavelet transform, output in pixel space.
%
% Default usage :
%
%   [f_wav, f_scal] = s2let_transform_analysis_mw(f, <options>)
%
% f is the input field -- MW sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'N'               = { Azimuthal/directional band-limit; N > 1 (default=L) }
%  'Spin'               = { Spin; (default=0) }
%  'J_min'           = { Minimum wavelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true       [full resolution wavelets] }
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               wavelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the wavelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the wavelets
%                       should be lowered from (default = 0)]
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
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
p.addParamValue('N', Lguessed, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(f, varargin{:});
args = p.Results;

if strcmp(args.Sampling, 'MWSS')
    f_vec = s2let_mwss_arr2vec(f);
else
    f_vec = s2let_mw_arr2vec(f);
end

if(all(isreal(f_vec)))
  f_vec = complex(f_vec,0);
end

[f_wav_vec, f_scal_vec] = s2let_transform_analysis_mw_mex(f_vec, args.B, args.L, args.J_min, ...
                                                          args.N, args.Spin, ...
                                                          args.Reality, args.Upsample, ...
                                                          args.SpinLowered, args.SpinLoweredFrom, ...
                                                          args.Sampling);

if strcmp(args.Sampling, 'MWSS')
    f_scal = s2let_mwss_vec2arr(f_scal_vec);

    J = s2let_jmax(args.L, args.B);
    f_wav = cell(J+1-args.J_min, args.N);
    offset = 0;
    for j = args.J_min:J
      for en = 1:args.N
        if args.Upsample
          band_limit = args.L;
        else
          band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
        end
        temp = zeros(band_limit+1, 2*band_limit);
        for t = 1:band_limit+1
            for p = 1:2*band_limit
                ind = offset + (t-1) * 2 * band_limit + p;
                temp(t,p) = f_wav_vec(1,ind);
            end
        end
        f_wav{j+1-args.J_min, en} = temp;
        offset = offset + (band_limit+1) * 2*band_limit;
      end
    end
else
    f_scal = s2let_mw_vec2arr(f_scal_vec);

    J = s2let_jmax(args.L, args.B);
    f_wav = cell(J+1-args.J_min, args.N);
    offset = 0;
    for j = args.J_min:J
      for en = 1:args.N
        if args.Upsample
          band_limit = args.L;
        else
          band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
        end
        temp = zeros(band_limit, 2*band_limit-1);
        for t = 1:band_limit
            for p = 1:2*band_limit-1
                ind = offset + (t-1) * ( 2 * band_limit - 1) + p;
                temp(t,p) = f_wav_vec(1,ind);
            end
        end
        f_wav{j+1-args.J_min, en} = temp;
        offset = offset + band_limit * (2*band_limit-1);
      end
    end
end
