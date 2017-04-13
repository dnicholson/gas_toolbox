% =========================================================================
% Fas_Fd: calculate diffusive air-sea flux for parameterizations without an
% explicit bubble flux
%
% [Fd, k] = fas_Fd(C,u10,S,T,slp,gas,param,rh)
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% [Fd, k] = fas_Fd(C,u10,S,T,slp,gas,param,rh)
% [Fd, k] = fas_Fd(0.01410,5,35,10,1,'Ar','W14',0.9)
%    > Fd = -4.1704e-09
%    > k = 1.6143e-05
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% C:        gas concentration (mol m-3)
% u10:      10 m wind speed (m/s)
% S:        Sea surface salinity (PSS)
% T:        Sea surface temperature (deg C)
% slp:      sea level pressure (atm)
% gas:      name of gas, formatted as a string, e.g. 'Ne'
% param     abbreviation for parameterization:
%           W14  = Wanninkhof 2014
%           W92a = Wanninkhof 1992 - averaged winds
%           W92b = Wanninkhof 1992 - instantaneous or steady winds
%           Sw07 = Sweeney et al. 2007
%           Ho06 = Ho et al. 2006
%           Ng00 = Nightingale et al. 2000
%           LM86 = Liss and Merlivat 1986
% rh:       relative humidity expressed as the fraction of saturation 
%           (0.5 = 50% RH).
%           rh is an optional but recommended argument. If not provided, it
%           will be set to 1 (100% RH) within the function.
%
%       Code    Gas name        Reference
%       ----   ----------       -----------
%       He      Helium          Weiss 1971
%       Ne      Neon            Hamme and Emerson 2004
%       Ar      Argon           Hamme and Emerson 2004
%       Kr      Krypton         Weiss and Keiser 1978
%       Xe      Xenon           Wood and Caputi 1966
%       N2      Nitrogen        Hamme and Emerson 2004   
%       O2      Oxygen          Garcia and Gordon 1992   
%
% Fd = -k.*(C-slpc.*Geq)
%
% Explanation of slpc:
%      slpc = (observed dry air pressure)/(reference dry air pressure)
% slpc is a pressure correction factor to convert from reference to
% observed conditions. Equilibrium gas concentration in gasmoleq is
% referenced to 1 atm total air pressure, including saturated water vapor
% (RH=1), but observed sea level pressure is usually different from 1 atm,
% and humidity in the marine boundary layer can be less than saturation.
% Thus, the observed sea level pressure of each gas will usually be
% different from the reference.
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% Fd        Diffusive flux                                (mol m-2 s-1)
% k:        Diffusive gas transfer velocity               (m s-1)
%
% -------------------------------------------------------------------------
% AUTHOR:
% -------------------------------------------------------------------------
% Author: Cara Manning cmanning@whoi.edu 
% Version: 12 April 2017
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2017 David Nicholson and Cara Manning 
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
% =========================================================================

function [Fd, k] = fas_Fd(C,u10,S,T,slp,gas,param,rh)

% if humidity is not provided, set to 1 for all values
if nargin == 8
    if mean(rh) < 0 || mean(rh) > 1
        error('Relative humidity must be 0 <= rh <= 1');
    end
else
    rh =1.*ones(size(C));
end

% slpc = (observed dry air pressure)/(reference dry air pressure)
% see Description section in header
ph2oveq = vpress(S,T);
ph2ov = rh.*ph2oveq;
slpc = (slp-ph2ov)./(1-ph2oveq);

[~,Sc] = gasmoldiff(S,T,gas);
Geq = gasmoleq(S,T,gas);
    
switch upper(param)
    case {'W14','W92A','W92B','SW07','HO06','NG00','LM86'}
        k = kgas(u10,Sc,param);
    otherwise
        error('Only W14, W92a, W92b, Sw07, Ho06, Ng00 and LM86 are supported.');
end

Fd = -k.*(C-slpc.*Geq);

end