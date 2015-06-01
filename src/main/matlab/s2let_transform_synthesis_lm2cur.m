function [flm_cur_syn, flm_scal_syn] = s2let_transform_synthesis_lm2cur(f_cur, f_scal, flm_gen, varargin)

% s2let_transform_synthesis_lm2cur
% Compute curvelet transform, input harmonic space output
% in pixel space.
%
% Default usage :
%
% [flm_cur_syn, flm_scal_syn] = s2let_transform_synthesis_lm2cur(f_cur, f_scal, <options>)
% flm_rec = s2let_transform_synthesis_lm2cur(f_cur, f_scal, <options>)
%
% f_cur contains the input curvelet contributions -- MW sampling,
% f_scal contains the input scaling contributions -- MW sampling,
% flm_gen : included here to compute Lguess from the input for "synthesis_lmn2lm",
% flm_rec is the output field -- harmonic space,
%
% Option :
%  'Reality'         = { false [do not assume corresponding signal f real (default)],
%                        true  [assume f real (improves performance)] }
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'N'               = { Azimuthal/directional band-limit; N > 1 (default=L) }
%  'Spin'               = { Spin; (default=0) }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true       [full resolution curvelets] }
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%  'J_min'           = { Minimum curvelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
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


%%%%%%%%%%%%
%{
len = size(f_cur);
temp = f_cur{len};
sz = size(temp);
if sz(1) == 2*sz(2)-1 || sz(2) == 2*sz(1)-1
    Lguessed = min([sz(1) sz(2)]);
else
    Lguessed = min([sz(1) sz(2)])-1;
end
%}
%%%%%%%%%%%%



p = inputParser;
p.addRequired('f_cur');
p.addRequired('f_scal', @isnumeric);
p.addRequired('flm_gen', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L',  @isnumeric);                             %Lguessed
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('N', @isnumeric);                              %Lguessed
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(f_cur, f_scal, flm_gen, varargin{:});
args = p.Results;


B= args.B;
L= args.L;
N= args.N;
Spin= args.Spin;
J_min= args.J_min;
J = s2let_jmax(args.L, args.B);


% ---------------
% Tile curvelets:
% ---------------
%***** step 1a) call curvelet- and scaling-function- generating functions
% disp('Tile curvelets in harmonic space (cur_lm, scal_l)')
[psi_lm phi_l] = s2let_curvelet_tiling(B, L, N, Spin, J_min);
for j = J_min:J,
cur_lm{j-J_min+1} = psi_lm(:,j+1);
end
%***** step 1b) compute the scaling coefficients (no j-dependence except on J_min)
scal_l = zeros(L^2,1);
for l = 0:L-1,
scal_l(l^2+l+1,1) = phi_l(l+1);
end



% -----------------
% Signal synthesis:
% -----------------

%***** step 1a) Transform to lmn space
%if strcmp(args.Sampling, 'MW')

% Scaling functions:
f_scal_lm_syn = ssht_forward(f_scal, L, 'spin', Spin);

% Curvelets:
for j = J_min:J,
f_cur_lmn_syn{j-J_min+1} = so3_forward(f_cur{j-J_min+1} , L, N);
end

%end

%***** step 1b) Reconstruct the signals in harmonic space
%***** step 1c) Do Wigner transform (lmn2lm -  Call matlab function synthesis_lmn2lm)
[flm_cur_syn, flm_scal_syn] = s2let_transform_synthesis_lmn2lm(f_cur_lmn_syn, f_scal_lm_syn, ...
                                                                cur_lm, scal_l, flm_gen, ...
                                                                'B', B, 'L', L, 'J_min', J_min, ...
                                                                'N', N, 'Upsample', false, 'Spin', 0);


%disp('Sum: flm_cur_syn+flm_scal_syn ');
%disp('then ');
%disp('Compute the re-constructed function via ssht_inverse ');
%flm_rec = flm_scal_syn+flm_cur_syn;


%%%%%%%%%%%%
%{
if  strcmp(args.Sampling, 'MWSS')
    f_scal_vec = s2let_mwss_arr2vec(f_scal);
else
    f_scal_vec = s2let_mw_arr2vec(f_scal);
end
if(all(isreal(f_scal_vec)))
  f_scal_vec = complex(f_scal_vec,0);
end
J = s2let_jmax(args.L, args.B);

f_cur_vec = [];

offset = 0;
if  strcmp(args.Sampling, 'MWSS')
    for j = args.J_min:J
      for en = 1:args.N
        if args.Upsample
            band_limit = args.L;
        else
            band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
        end
        temp = f_cur{j+1-args.J_min, en};
        for t = 1:band_limit+1
            for p = 1:2*band_limit
               ind = offset + (t-1) * 2 * band_limit + p;
                f_cur_vec = [f_cur_vec temp(t,p)];
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
          temp = f_cur{j+1-args.J_min, en};
          for t = 1:band_limit
              for p = 1:2*band_limit-1
                ind = offset + (t-1) * ( 2 * band_limit - 1) + p;
                f_cur_vec = [f_cur_vec temp(t,p)];
              end
          end
          offset = offset + band_limit * (2 * band_limit - 1);
      end
    end
end


if(all(isreal(f_cur_vec)))
  f_cur_vec = complex(f_cur_vec,0);
end

    
%%%%%
flm = s2let_transform_synthesis_lm2cur_mex(f_cur_vec, f_scal_vec, args.B, args.L, args.J_min, ...
                                           args.N, args.Spin, args.Reality, args.Upsample, ...
                                           args.SpinLowered, args.SpinLoweredFrom, ...
                                           args.Sampling);

end
%}
%%%%%%%%%%%%%
