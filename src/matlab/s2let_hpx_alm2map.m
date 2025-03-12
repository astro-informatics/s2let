function f = s2let_hpx_alm2map(flm, nside, varargin)

% s2let_axisym_synthesis 
% Compute axisymmetric wavelet transform, output as HEALPIX maps.
%
% Default usage :
%
%   f = s2let_hpx_alm2map(flm, nside, <options>)
%
% flm is the input spherical harmonic decomposition,
% nside is the HEALPIX resolution for the output map,
% f is the corresponding output map.
%
% Option :
%  'L'               = { Harmonic band-limit; L > 1 (default=2*nside) }
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details


sz = size(flm);
Lguessed = sqrt(sz(1));

p = inputParser;
p.addRequired('flm',@isnumeric); 
p.addRequired('nside', @isnumeric); 
p.addParamValue('L', Lguessed, @isnumeric);   
p.parse(flm, nside, varargin{:});
args = p.Results;

f = s2let_hpx_alm2map_mex(flm, nside, args.L);
 
end