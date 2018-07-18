function [flmQ, flmU] = s2let_hpx_map2alm_spin(fQ, fU, varargin)
% s2let_hpx_axisym_analysis_spin
%
% Default usage :
%
%   flmQ, flmU = s2let_hpx_map2alm_spin(fQ, fU, <options>)
%
% f is the input field -- HEALPIX sampling,
% flm is the output spherical harmonic decomposition.
%
% Option :
%  'nside'           = { HEALPIX resolution; (default=guessed)}
%  'Spin'           = { Spin; (default=2)}
%  'L'               = { Harmonic band-limit; L > 1 (default=2*nside) }

% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(fQ);
nsideguessed = sqrt(max(sz)/12);
Lguessed = 2*nsideguessed;

p = inputParser;
p.addRequired('fQ', @isnumeric); 
p.addRequired('fU', @isnumeric); 
p.addParamValue('nside', nsideguessed, @isnumeric);  
p.addParamValue('Spin', 2, @isnumeric);    
p.addParamValue('L', Lguessed, @isnumeric);   
p.parse(fQ, fU, varargin{:});
args = p.Results;

[flmQ, flmU] = s2let_hpx_map2alm_spin_mex(fQ, fU, args.nside, args.L, args.Spin);

end
