function flm = s2let_transform_synthesis_lm2wav(f_wav, f_scal, varargin)

% s2let_transform_synthesis_lm2wav
% Compute spin directional wavelet transform, input harmonic space output
% in pixel space.
%
% Default usage :
%
%   flm = s2let_transform_synthesis_lm2wav(f_wav, f_scal, <options>)
%
% f_wav contains the input wavelet contributions -- MW sampling,
% f_scal contains the input scaling contributions -- MW sampling,
% flm is the output field -- harmonic space,
%
% Option :
%  'Reality'         = { false [do not assume corresponding signal f real (default)],
%                        true  [assume f real (improves performance)] }
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'N'               = { Azimuthal/directional band-limit; N > 1 (default=L) }
%  'Spin'               = { Spin; (default=0) }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true       [full resolution wavelets] }
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%  'J_min'           = { Minimum wavelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'OriginalSpin' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the wavelets
%                       should be lowered from (default = 0)]
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

len = size(f_wav);
temp = f_wav{len};
sz = size(temp);
if sz(1) == 2*sz(2)-1 || sz(2) == 2*sz(1)-1
    Lguessed = min([sz(1) sz(2)]);
else
    Lguessed = min([sz(1) sz(2)])-1;
end

p = inputParser;
p.addRequired('f_wav');
p.addRequired('f_scal', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('N', Lguessed, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('OriginalSpin', 0, @isnumeric);
p.parse(f_wav, f_scal, varargin{:});
args = p.Results;

if  strcmp(args.Sampling, 'MWSS')
    f_scal_vec = s2let_mwss_arr2vec(f_scal);
else
    f_scal_vec = s2let_mw_arr2vec(f_scal);
end
if(all(isreal(f_scal_vec)))
  f_scal_vec = complex(f_scal_vec,0);
end
J = s2let_jmax(args.L, args.B);

f_wav_vec = [];

offset = 0;
if  strcmp(args.Sampling, 'MWSS')
    for j = args.J_min:J
      for en = 1:args.N
        if args.Upsample
            band_limit = args.L;
        else
            band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
        end
        temp = f_wav{j+1-args.J_min, en};
        for t = 1:band_limit+1
            for p = 1:2*band_limit
               ind = offset + (t-1) * 2 * band_limit + p;
                f_wav_vec = [f_wav_vec temp(t,p)];
            end
        end
        offset = offset + (band_limit+1) * 2 * band_limit;
      end
    end
else
    for j = args.J_min:J
      for en = 1:args.N
        if args.Upsample
            band_limit = args.L;
          else
            band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
          end
          temp = f_wav{j+1-args.J_min, en};
          for t = 1:band_limit
              for p = 1:2*band_limit-1
                ind = offset + (t-1) * ( 2 * band_limit - 1) + p;
                f_wav_vec = [f_wav_vec temp(t,p)];
              end
          end
          offset = offset + band_limit * (2 * band_limit - 1);
      end
    end
end


if(all(isreal(f_wav_vec)))
  f_wav_vec = complex(f_wav_vec,0);
end

flm = s2let_transform_synthesis_lm2wav_mex(f_wav_vec, f_scal_vec, args.B, args.L, args.J_min, ...
                                           args.N, args.Spin, args.Reality, args.Upsample, ...
                                           args.OriginalSpin, ...
                                           args.Sampling);

end
