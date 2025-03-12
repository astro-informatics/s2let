function map_sm = s2let_smoothing(map, L, fwhm, varargin)

p = inputParser;
p.addParamValue('Reality', false, @islogical);
p.addParamValue('Method', 'MW', @ischar);
p.parse(varargin{:});
args = p.Results;

map_lm = ssht_forward(map, L, 'Reality', true, 'Method', args.Method);
map_sm_lm = zeros(size(map_lm));
sigma = fwhm / 2.355;
    
ind = 1;
for el = 0:L-1
    for m = -el:el
        map_sm_lm(ind) = map_lm(ind) * exp(-el*(el+1) * sigma^2 / 2);
        ind = ind+1;
    end
end

map_sm = ssht_inverse(map_sm_lm, L, 'Reality', args.Reality, 'Method', args.Method);
