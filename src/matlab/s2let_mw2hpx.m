function f_hpx = s2let_mw2hpx(f, nside, varargin)

% s2let_mw2hpx 
% Converts MW map into Healpix map
%
% Default usage :
%
%   f_hpx = s2let_hpx2mw(f_mw, <options>)
%
% f is the input map -- MW sampling,
% f is the output map -- Healpix sampling.
%
% Option :
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed) }

% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
Lguessed = min([sz(1) sz(2)]);

p = inputParser;
p.addRequired('f', @isnumeric); 
p.addParamValue('L', Lguessed, @isnumeric);   
p.parse(f, varargin{:});
args = p.Results;

[flm] = ssht_forward(f, args.L, 'Reality', true);
f_hpx = s2let_hpx_alm2map(flm, nside, 'L', args.L);

end
