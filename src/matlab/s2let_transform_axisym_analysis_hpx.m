function [f_wav, f_scal] = s2let_transform_axisym_analysis_hpx(f, varargin)

% s2let_transform_axisym_analysis_hpx 
% Compute axisymmetric wavelet transform, output as HEALPIX maps.
%
% Default usage :
%
%   [f_wav, f_scal] = s2let_transform_axisym_analysis_hpx(f, <options>)
%
% f is the input field -- HEALPIX sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
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

sz = size(f);
nsideguessed = sqrt(max(sz)/12);
Lguessed = 2*nsideguessed;

p = inputParser;
p.addRequired('f', @isnumeric); 
p.addParamValue('nside', nsideguessed, @isnumeric);   
p.addParamValue('B', 2, @isnumeric);   
p.addParamValue('L', Lguessed, @isnumeric);   
p.addParamValue('J_min', 0, @isnumeric); 
p.parse(f, varargin{:});
args = p.Results;

[f_wav_vec, f_scal] = s2let_transform_axisym_analysis_hpx_mex(f, args.nside, args.B, args.L, args.J_min);

J = s2let_jmax(args.L, args.B);
npix = 12*args.nside*args.nside;
f_wav = cell(J+1-args.J_min, 1);
offset = 1;
for j = args.J_min:J
  f_wav{j+1-args.J_min} = f_wav_vec(offset:offset+npix-1);
  offset = offset + npix;
end
