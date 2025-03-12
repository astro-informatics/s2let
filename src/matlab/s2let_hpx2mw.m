function f_mw = s2let_hpx2mw(f, varargin)

% s2let_hpx2mw 
% Converts Healpix map into MW map
%
% Default usage :
%
%   f_mw = s2let_hpx2mw(f_hpx, <options>)
%
% f is the input map -- HEALPIX sampling,
% f is the output map -- MW sampling.
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

flm = s2let_hpx_map2alm(f, 'nside', args.nside, 'L', args.L);
%flm = s2let_hpx_map2alm_mex(f, args.nside, args.L);
f_mw = ssht_inverse(flm, args.L, 'Reality', true);

end
