function s2let_plot_sphere(wav, scal, B, L, N, J_min, varargin)
% s2let_plot_sphere
% Plot wavelet coefficients on multiple spheres. The function
% generates one plot of the scaling function contribution and
% a grid of plots for each orientation of each scale of the
% wavelet contributions.
%
% Default usage :
%
%   s2let_plot_sphere(wav, scal, B, L, N, J_min, <options>)
%
% wav is cell array with all the wavelet coefficients.
% its first index is the wavelet scale j, the second
% index is the orientation g, and each element is a
% function on the sphere in MW sampling.
% scal is the corresponding scaling function contribution
% (i.e. just a single function on the sphere).
% B is the wavelet parameter.
% L is the angular band-limit.
% N is the orientational band-limit.
% J_min is the first wavelet scale in wav.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'Upsample'      = { false   [multiresolution algorithm (default)],
%                      true  [full resolution wavelets] },
%  'Function'        = { 'real' [plot the real part of the input functions (default)],
%                        'imag' [plot the imaginary part of the input functions],
%                        'abs'  [plot the absolute value of the input functions] }

% TODO: make more parameters optional and guess them from the sizes of
% wav and scal

% Parse arguments.
p = inputParser;
p.addRequired('wav');
p.addRequired('scal', @isnumeric);
p.addRequired('B', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('N', @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Function', 'real', @ischar)

p.parse(wav, scal, B, L, N, J_min, varargin{:});

args = p.Results;

J = s2let_jmax(L,B);

switch args.Function
    case {'real', 'imag', 'abs'}
        f = str2func(args.Function);
    otherwise
        error('Invalid Function parameter. Must be one of "real", "imag", "abs".')
end

figure(1);

if args.Upsample
    bl = L;
else
    bl = min([s2let_bandlimit(J_min-1, J_min, B, L), L]);
end
ssht_plot_sphere(f(scal), bl)

figure(2);

iplot = 1;
for j = J_min:J
    for n = 1:2*N-1
        subplot(J-J_min+1,2*N-1,iplot);
        if args.Upsample
            bl = L;
        else
            bl = min([s2let_bandlimit(j, J_min, B, L), L]);
        end
        ssht_plot_sphere(f(wav{j-J_min+1,n}), bl)
        iplot = iplot + 1;
    end
end
