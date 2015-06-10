function f_rec = s2let_transform_curvelet_mw(f_init, varargin)

% s2let_transform_analysis_mw
% Perform curvelet transform in harmonic space and then reconstruction (output in pixel space).
%  
% Default usage :
%
%   [f_rec] = s2let_transform_analysis_mw(f_init, <options>)
%
% f_init is the input field -- MW sampling,
% f_rec is the output reconstructed field. 
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
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

sz = size(f_init);
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
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(f_init, varargin{:});
args = p.Results;
N = args.L; 


[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_init,...
                                                         'B',args.B, 'L', args.L, 'J_min', args.J_min, ...
                                                         'Spin', args.Spin,'Reality', args.Reality,...
                                                         'Upsample', args.Upsample, ...
                                                         'SpinLowered', args.SpinLowered, ...
                                                         'SpinLoweredFrom',  args.SpinLoweredFrom, ...
                                                         'Sampling', args.Sampling);

f_rec = s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                  'B',args.B, 'L', args.L, 'J_min', args.J_min, ...
                                                  'Spin', args.Spin,'Reality', args.Reality,...
                                                  'Upsample', args.Upsample, ...
                                                  'SpinLowered', args.SpinLowered, ...
                                                  'SpinLoweredFrom',  args.SpinLoweredFrom, ...
                                                  'Sampling', args.Sampling);

                                                      
% ====================
% Below shows the equivalent version of curvelet transform: 
% ====================
%{
% ---------------
% Transform the signals in harmonic space (ssht_forward):
% ---------------
% disp('Construct the corresponding flm of the signal via ssht_forward'); 
flm_init= ssht_forward(f_init, args.L, 'Method', args.Sampling);

% ---------------
% Tile curvelets:
% ---------------
% disp('curvelet_tiling: Tile curvelets in harmonic space (cur_lm, scal_l)')
[cur_lm scal_l] = s2let_curvelet_tiling(args.B, args.L, args.J_min, ...
                                        'Spin', args.Spin, 'SpinLowered', args.SpinLowered,...
                                        'SpinLoweredFrom',args.SpinLoweredFrom);
% where 
% cur_lm, @iscell, contains the curvelet kernels in harmonic space,
% scal_l, @isnumeric, contains the scaling contributions in harmonic space


% -----------------
% Signal analysis: 
% -----------------
% disp('analysis_lm2lmn: from the harmonic to the Wigner space ')
[f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_init, cur_lm, scal_l, ...
                                                         'B', args.B, 'L', args.L, 'J_min', args.J_min,...
                                                         'Spin', args.Spin,  ...
                                                         'Reality', args.Reality, 'Upsample', args.Upsample,...
                                                         'SpinLowered',args.SpinLowered,  'SpinLoweredFrom', args.SpinLoweredFrom, ...
                                                         'Sampling',  args.Sampling);

                                                     
J =s2let_jmax(args.L, args.B);
% disp('so3_inverse: f_cur_lmn to f_cur ; then so3_forward to get back f_cur_lmn')
for j = args.J_min:J,
f_cur{j-args.J_min+1} = so3_inverse(f_cur_lmn{j-args.J_min+1}, args.L, N, 'Sampling', args.Sampling,'Reality', args.Reality);
f_cur_lmn{j-args.J_min+1} = so3_forward(f_cur{j-args.J_min+1}, args.L, N, 'Sampling', args.Sampling,'Reality', args.Reality);
end  
% scaling function
f_scal= ssht_inverse(f_scal_lm, args.L,'Method', args.Sampling);
f_scal_lm = ssht_forward(f_scal, args.L,'Method', args.Sampling);                                                     

% -----------------
% Signal synthesis: (pixel to harmonic space) 
% -----------------
%  disp('synthesis_lmn2lm: reconstruct the signal in hamonic space') 
flm_rec  = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l, ...
                                                    'B', args.B, 'L', args.L, 'J_min', args.J_min,...
                                                    'Spin', args.Spin,  ...
                                                    'Reality', args.Reality, 'Upsample', args.Upsample,...
                                                    'SpinLowered',args.SpinLowered,  'SpinLoweredFrom', args.SpinLoweredFrom, ...
                                                    'Sampling',  args.Sampling);                                                     

% disp('ssht_inverse : Reconstruct the siganls form their harmonic coefficients');
f_rec = ssht_inverse(flm_rec, args.L, 'Method', args.Sampling); 
%}                                                
                                                     
                                           
                                                      
end