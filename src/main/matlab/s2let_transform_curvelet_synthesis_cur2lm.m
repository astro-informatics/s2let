function flm_rec = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal,  varargin)

% s2let_transform_synthesis_lm2cur
% Compute curvelet transform:  
% input in curvelet space (i.e. hamonic to Wigner space via SO3_forward) 
% output in harmonic space (i.e. Wigner to harmonic space via synthesis_lmn2lm) .
%
% Default usage :
%
% flm_rec = s2let_transform_curvelet_synthesis_lm2cur(f_cur, f_scal,  <options>)
%
% f_cur contains the input curvelet contributions -- MW sampling,
% f_scal contains the input scaling contributions -- MW sampling,
% flm_rec is the output field = flm_cur_syn+ flm_scal_syn
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'J_min'           = { Minimum curvelet scale to consider;
%  'Spin'               = { Spin; (default=0) }
%                        0 <= J_min < log_B(L) (default=0) }
%  'Reality'         = { false        [do not assume corresponding signal f real (default)],
%                        true         [assume f real (improves performance)] }
%  'Upsample'        = { false        [multiresolution algorithm (default)],
%                        true       [full resolution curvelets] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               curvelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the curvelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the curvelets
%                       should be lowered from (default = 0)]
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%
% -----------------------------------------------------------
% Log:
% -  constructed by Jennifer Y H Chan on 5th June 2015
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------



len = size(f_cur);
temp = f_cur{len};
sz = size(temp);
if sz(1) == 2*sz(2)-1 || sz(2) == 2*sz(1)-1
    Lguessed = min([sz(1) sz(2)]);
else
    Lguessed = min([sz(1) sz(2)])-1;
end


p = inputParser;
p.addRequired('f_cur');
p.addRequired('f_scal', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);                           
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.addParamValue('Sampling', 'MW', @ischar);
p.parse(f_cur, f_scal, varargin{:});
args = p.Results;


N= args.L;
Nj=N;
band_limit = args.L;
J = s2let_jmax(args.L, args.B);

% ---------------
% Tile curvelets:
% ---------------
% Call curvelet- and scaling-function- generating functions
[cur_lm scal_l] = s2let_curvelet_tiling(args.B, args.L, args.J_min,...
                                        'Spin', args.Spin,...
                                        'SpinLowered', args.SpinLowered, ...
                                        'SpinLoweredFrom', args.SpinLoweredFrom);



% -----------------
% Signal synthesis: (Transform to lmn space, then reconstruct the signal in harmonic space)
% -----------------
% Scaling functions:
f_scal_lm_syn = ssht_forward(f_scal, args.L, 'Spin', args.Spin, ...
                            'Method', args.Sampling, 'Reality', args.Reality);
% Curvelet contribution:
for j = args.J_min:J,
 if (args.Upsample ==0)  %false => multi-resolution 
     band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
     Nj = min(N, band_limit);
     Nj = Nj+ mod((Nj+N),2) ;  
 end
f_cur_lmn_syn{j-args.J_min+1} = so3_forward(f_cur{j-args.J_min+1} ,band_limit, Nj, ...
                                            'Reality', args.Reality, 'Sampling', args.Sampling);
end


% Reconstruct the signals in harmonic space
% Perform Wigner transform (lmn2lm -  Call matlab function synthesis_lmn2lm)
flm_rec = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn_syn, f_scal_lm_syn, cur_lm, scal_l,...
                                                    'B', args.B, 'L', args.L, ...
                                                    'Spin', args.Spin, ...
                                                    'J_min', args.J_min, ...
                                                    'Upsample', args.Upsample,...
                                                    'Reality', args.Reality,...
                                                    'SpinLowered', args.SpinLowered, ...
                                                    'SpinLoweredFrom', args.SpinLoweredFrom,...
                                                    'Sampling', args.Sampling );

                                      
end