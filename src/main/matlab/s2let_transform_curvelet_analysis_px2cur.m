function [f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_init, varargin)

% s2let_transform_analysis_lm2cur
% Compute spin directional curvelet transform, input in harmonic space,
% output in pixel space.
%
% Default usage :
%
%   [f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_init, <options>)
%
% f_init is the input field in pixel space,
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
% S2LET package to perform curvelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f_init);
if sz(1) == 2*sz(2)-1 || sz(2) == 2*sz(1)-1
    Lguessed = min([sz(1) sz(2)]);
else
    Lguessed = min([sz(1) sz(2)])-1;
end

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
p.parse(f_init, varargin{:});
args = p.Results;

% For curvelets, azimuthal/directional band-limit N always equals to L           
N = args.L ; 
J = s2let_jmax(args.L, args.B);

% ---------------
% Construct signals in harmonic space:
% ---------------
flm_init= ssht_forward(f_init, args.L,  'Method', args.Sampling, 'Spin', args.Spin, 'Reality', args.Reality);

% ---------------
% Signal analysis (from harmonic to curvelet space):
% ---------------
[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_init,  ...
                                                         'B',args.B, 'L', args.L, 'J_min', args.J_min, ...
                                                         'Spin', args.Spin,'Reality', args.Reality,...
                                                         'Upsample', args.Upsample, ...
                                                         'SpinLowered', args.SpinLowered, ...
                                                         'SpinLoweredFrom',  args.SpinLoweredFrom, ...
                                                         'Sampling', args.Sampling);
                                                     
% Clear arrary memory:
flm_init = 0;

end