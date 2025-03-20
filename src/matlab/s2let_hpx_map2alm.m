function flm = s2let_hpx_map2alm(f, varargin)

% s2let_hpx_axisym_analysis 
% Compute axisymmetric wavelet transform, output as HEALPIX maps.
%
% Default usage :
%
%   flm = s2let_hpx_map2alm(f, <options>)
%
% f is the input field -- HEALPIX sampling,
% flm is the output spherical harmonic decomposition.
%
% Option :
%  'nside'           = { HEALPIX resolution; (default=guessed)}
%  'L'               = { Harmonic band-limit; L > 1 (default=2*nside) }

% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
nsideguessed = sqrt(max(sz)/12);
Lguessed = 2*nsideguessed;

p = inputParser;
p.addRequired('f', @isnumeric); 
p.addParamValue('nside', nsideguessed, @isnumeric);   
p.addParamValue('L', Lguessed, @isnumeric);   
p.parse(f, varargin{:});
args = p.Results;

flm = s2let_hpx_map2alm_mex(f, args.nside, args.L);

end