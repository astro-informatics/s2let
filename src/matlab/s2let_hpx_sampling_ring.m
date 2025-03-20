function [thetas, phis] = s2let_hpx_sampling_ring(nside)

% s2let_hpx_sampling_ring 
% Compute the Healpix sampling nodes
% Default usage :
%
%   [thetas, phis] = s2let_hpx_sampling_ring(nside)
%
% nside is the input Healpix resolution,
% thetas contains the colatitudes of the nodes,
% phis contains the longitudes of the nodes.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

%     npix = 12 * nside^2;
%     thetas = zeros(npix,1);
%     phis = zeros(npix,1);
%     for i = 1:npix
%        [t,p] = healpix_pix2ang_ring(i, nside);
%        thetas(i) = t;
%        phis(i) = p;
%     end

    npix = 12 * nside^2;
    ipix = 1:npix;
    nl2 = 2 * nside;
    ncap = 2  * nside * (nside - 1);  % points in each polar cap, =0 for nside =1

    thetas = zeros(npix,1);
    phis = zeros(npix,1);
    
    ind = (ipix <= ncap); % North polar cap
        hip   = ipix(ind) .* 0.5;
        fihip = fix(hip);
        iring = fix( sqrt( hip - sqrt(fihip) ) ) + 1 ;% counted from North pole
        iphi  = ipix(ind) - 2*iring.*(iring - 1) ;

        thetas(ind) = acos( 1.0 - iring.^2.0 ./ (3.0*nside^2.0) );
        phis(ind)   = (iphi - 0.5) .* pi ./ (2.0.*iring);

    ind = (ipix <= nl2*(5*nside+1) & ipix > ncap); % Equatorial region

       ip    = ipix(ind) - ncap - 1;
       nl4   = 4*nside;
       iring = floor( ip / nl4 ) + nside; % counted from North pole
       iphi  = mod(ip,nl4) + 1;

       fodd  = 0.5 * (1 + mod(iring+nside,2)) ; % 1 if iring+nside is odd, 1/2 otherwise
       thetas(ind) = acos( (nl2 - iring) ./ (1.5*nside) );
       phis(ind)   = (iphi - fodd) * pi /(2.0*nside);

    ind = (ipix > nl2*(5*nside+1) & ipix > ncap); % South polar cap

       ip    = npix - ipix(ind) + 1;
       hip   = ip*0.5;
       fihip = round(hip);
       iring = floor( sqrt( hip - sqrt(fihip) ) ) + 1 ; % counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring.*(iring-1));

       thetas(ind) = acos( -1.0 + iring.^2 ./ (3.0*nside^2) );
       phis(ind)   = (iphi - 0.5) * pi ./ (2.0*iring);

end

function [theta, phi] = healpix_pix2ang_ring(ipix, nside)

    npix = 12 * nside^2;
    nl2 = 2 * nside;
    ncap = 2  * nside * (nside - 1);  % points in each polar cap, =0 for nside =1

    if (ipix <= ncap)
        hip   = ipix * 0.5;
        fihip = fix(hip);
        iring = fix( sqrt( hip - sqrt(fihip) ) ) + 1 ;% counted from North pole
        iphi  = ipix - 2*iring*(iring - 1) ;

        theta = acos( 1.0 - iring^2.0 / (3.0*nside^2.0) );
        phi   = (iphi - 0.5) * pi/(2.0*iring);

    elseif (ipix <= nl2*(5*nside+1)) % Equatorial region ------

       ip    = ipix - ncap - 1;
       nl4   = 4*nside;
       iring = floor( ip / nl4 ) + nside; % counted from North pole
       iphi  = mod(ip,nl4) + 1;

       fodd  = 0.5 * (1 + mod(iring+nside,2)) ; % 1 if iring+nside is odd, 1/2 otherwise
       theta = acos( (nl2 - iring) / (1.5*nside) );
       phi   = (iphi - fodd) * pi /(2.0*nside);

    else % South Polar cap -----------------------------------

       ip    = npix - ipix + 1;
       hip   = ip*0.5;
       fihip = round(hip);
       iring = floor( sqrt( hip - sqrt(fihip) ) ) + 1 ; % counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

       theta = acos( -1.0 + iring^2 / (3.0*nside^2) );
       phi   = (iphi - 0.5) * pi/(2.0*iring);

    end

end