function [f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_init, varargin)

% s2let_transform_analysis_lm2cur
% Compute spin directional curvelet transform, input in harmonic space,
% output in pixel space.
%
% Default usage :
%
%   [f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_init, <options>)
%
% flm_init is the input field in harmonic space,
% f_cur contains the output curvelet contributions,
% f_scal contains the output scaling contributions,
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

sz = length(flm_init(:));
Lguessed = sqrt(sz);

p = inputParser;
p.addRequired('flm_init', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.addParamValue('Sampling', 'MW', @ischar);
p.parse(flm_init, varargin{:});
args = p.Results;

% For curvelets, azimuthal/directional band-limit N always equals to L           
N = args.L ; 
% Nj=N;
% band_limit = args.L;
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
% Signal analysis:
% -----------------
% Decompose the signals using curvelets and the scaling functions
% And perform Wigner transform (lm2lmn -  Call matlab function s2let_transform_analysis_lm2lmn)
[f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_init, cur_lm, scal_l,...
                                                                'B',args.B, 'L', args.L, 'J_min', args.J_min, ...
                                                                'Spin', args.Spin,'Reality', args.Reality,...
                                                                'Upsample', args.Upsample, ...
                                                                'SpinLowered', args.SpinLowered, ...
                                                                'SpinLoweredFrom',  args.SpinLoweredFrom, ...
                                                                'Sampling', args.Sampling);

% -----------------                                                     
% Transform to pixel space:
% -----------------
% Scaling functions:
f_scal = ssht_inverse(f_scal_lm, args.L, 'Spin', args.Spin, ...
                      'Method', args.Sampling, 'Reality', args.Reality);                  

 %{
 if (args.Upsample == 0)  %false => multi-resolution 
     band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
     Nj = band_limit; 
    % Nj = min(N, band_limit);
    % Nj = Nj+ mod((Nj+N),2) ;  %ensure N and Nj are both even and both odd
 end
%}  
% Curvelets:
for j = args.J_min:J,  
 f_cur{j-args.J_min+1} = so3_inverse(f_cur_lmn{j-args.J_min+1}, args.L, N , ...
                        'Sampling', args.Sampling,'Reality', args.Reality) ;
end
% band_limit, Nj, ...


flm_int = 0;

end