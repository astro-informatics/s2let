function [f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm, varargin)

% s2let_transform_analysis_lm2cur
% Compute spin directional curvelet transform, input in harmonic space,
% output in pixel space.
%
% Default usage :
%
%   [f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm, <options>)
%
% flm is the input field in harmonic space,
% f_cur contains the output curvelet contributions,
% f_scal contains the output scaling contributions,
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'N'               = { Azimuthal/directional band-limit; N > 1 (default=L) }
%  'Spin'               = { Spin; (default=0) }
%  'J_min'           = { Minimum curvelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true       [full resolution curvelets] }
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%  'Reality'         = { false        [do not assume corresponding signal f real (default)],
%                        true         [assume f real (improves performance)] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               curvelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the curvelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the curvelets
%                       should be lowered from (default = 0)]
%
% S2LET package to perform curvelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = length(flm(:));
Lguessed = sqrt(sz);

p = inputParser;
p.addRequired('flm', @isnumeric);
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
p.parse(flm, varargin{:});
args = p.Results;

flm_vec = flm(:);

if(all(isreal(flm_vec)))
  flm_vec = complex(flm_vec,0);
end

[f_cur_vec, f_scal_vec] = s2let_transform_analysis_lm2cur_mex(flm_vec, args.B, args.L, args.J_min, ...
                                                              args.N, args.Spin, ...
                                                              args.Reality, args.Upsample, ...
                                                              args.SpinLowered, args.SpinLoweredFrom, ...
                                                              args.Sampling);

if strcmp(args.Sampling, 'MWSS')
    f_scal = s2let_mwss_vec2arr(f_scal_vec);

    J = s2let_jmax(args.L, args.B);
    f_cur = cell(J+1-args.J_min, args.N);
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
                temp(t,p) = f_cur_vec(1,ind);
            end
        end
        f_cur{j+1-args.J_min, en} = temp;
        offset = offset + (band_limit+1) * 2*band_limit;
      end
    end
else
    f_scal = s2let_mw_vec2arr(f_scal_vec);

    J = s2let_jmax(args.L, args.B);
    f_cur = cell(J+1-args.J_min, args.N);
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
                temp(t,p) = f_cur_vec(1,ind);
            end
        end
        f_cur{j+1-args.J_min, en} = temp;
        offset = offset + band_limit * (2*band_limit-1);
      end
    end
end
