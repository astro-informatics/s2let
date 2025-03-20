function f = s2let_transform_axisym_synthesis_hpx(f_wav, f_scal, varargin)

% s2let_axisym_synthesis 
% Compute axisymmetric wavelet transform, input/outputs as HEALPIX maps.
%
% Default usage :
%
%   f = s2let_transform_axisym_synthesis_hpx(f_wav, f_scal, <options>)
%
% f_wav contains the input wavelet contributions -- HEALPIX sampling,
% f_scal contains the input scaling contributions -- HEALPIX sampling,
% f is the output field -- HEALPIX sampling,
%
% Option :
%  'nside'           = { HEALPIX resolution; (default=guessed)}
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=2*nside) }
%  'J_min'           = { Minimum wavelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details


len = size(f_wav);
temp = f_wav{len};
sz = size(temp);
nsideguessed = sqrt(max(sz)/12);
Lguessed = 2*nsideguessed;

p = inputParser;
p.addRequired('f_wav'); 
p.addRequired('f_scal', @isnumeric); 
p.addParamValue('nside', nsideguessed, @isnumeric); 
p.addParamValue('B', 2, @isnumeric);   
p.addParamValue('L', 2*nsideguessed, @isnumeric);   
p.addParamValue('J_min', 0, @isnumeric); 
p.parse(f_wav, f_scal, varargin{:});
args = p.Results;


J = s2let_jmax(args.L, args.B);
f_wav_vec = [];
for j = args.J_min:J
    temp = f_wav{j+1-args.J_min};
    f_wav_vec = [f_wav_vec temp];
end

f = s2let_transform_axisym_synthesis_hpx_mex(f_wav_vec, f_scal, args.nside, args.B, args.L, args.J_min);
 
end
