function f_rec = s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal,  varargin)

% s2let_transform_curvelet_synthesis_cur2px
% Compute (spin) curvelet transform,
% input curvelet space
% output in pixel space.
%
% Default usage :
%
% f_rec = s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal,  flm_init, <options>)
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
%
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

len = length(f_cur);   % i.e. f_cur(Nj-1 Nj, 2L-1)
temp = f_cur{len};
sz = size(temp);
if sz(1) == 2*sz(2)-1 || sz(1) == sz(3)
    Lguessed = sz(2); 
else
    Lguessed = (sz(3)+1)/2 ;
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

J = s2let_jmax(args.L, args.B);

% -----------------
% Signal synthesis:
% -----------------
% Reconstruct the signals in harmonic space:
flm_rec = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal,  ...
                                                    'B', args.B, 'L', args.L, ...
                                                    'Spin', args.Spin, ...
                                                    'J_min', args.J_min, ...
                                                    'Upsample', args.Upsample,...
                                                    'Reality', args.Reality,...
                                                    'SpinLowered', args.SpinLowered, ...
                                                    'SpinLoweredFrom', args.SpinLoweredFrom,...
                                                    'Sampling', args.Sampling );
                                      
% Reconstruct the signals in pxiel space:   
f_rec = ssht_inverse(flm_rec, args.L, ...
                    'Spin', args.Spin, ... 
                    'Method', 'MW', ...
                    'Reality', args.Reality);      
                                      
% Clear array memory:                                    
flm_rec = 0; 

end