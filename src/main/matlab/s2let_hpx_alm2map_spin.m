function [fQ, fU] = s2let_hpx_alm2map(flmQ, flmU, nside, varargin)

% s2let_axisym_synthesis_spin
% Compute axisymmetric wavelet transform, output as HEALPIX maps.
%
% Default usage :
%
%   fQ, fU = s2let_hpx_alm2map(flm, nside, <options>)
%
% flm is the input spherical harmonic decomposition,
% nside is the HEALPIX resolution for the output map,
% f is the corresponding output map.
%
% Option :
%  'L'               = { Harmonic band-limit; L > 1 (default=2*nside) }
%  'Spin'           = { Spin; (default=2)}
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details


sz = size(flmQ);
Lguessed = sqrt(sz(1));

p = inputParser;
p.addRequired('flmQ',@isnumeric); 
p.addRequired('flmU',@isnumeric); 
p.addRequired('nside', @isnumeric); 
p.addParamValue('Spin', 2, @isnumeric);    
p.addParamValue('L', Lguessed, @isnumeric);   
p.parse(flmQ, flmU, nside, varargin{:});
args = p.Results;

[fQ, fU] = s2let_hpx_alm2map_spin_mex(flmQ, flmU, nside, args.L, args.Spin);
 
end
