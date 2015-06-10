function s2let_spin_lowered_norm_factor = s2let_spin_lowered_normalization(el, varargin)

% s2let_spin_lowered_normalization(l, Spin)
% Computes the normalization factor for spin-lowered wavelets,
% which is sqrt((l+s)!/(l-s)!).
%
% Default usage :
%
%   s2let_spin_lowered_normalization(l, original_spin)
%
% l is the input harmonic index l.
% Spin is the spin number that the wavelet was lowered from.
%
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2015  
% -----------------------------------------------------------
% S2LET package to perform curvelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details


p = inputParser;
p.addRequired('el', @isnumeric);
p.addParamValue('original_spin', 0, @isnumeric);
p.parse(el, varargin{:});
args = p.Results;


factor = 1.0;
s = 0;
for s = (-abs(args.original_spin)+1):(abs(args.original_spin))
    factor = factor*(args.el+s);
end

if (args.original_spin > 0)
s2let_spin_lowered_norm_factor = sqrt(factor);
else
s2let_spin_lowered_norm_factor = sqrt(1.0/factor);
end



end

