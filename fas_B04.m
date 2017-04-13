% fas_BO4: Function to calculate air-sea fluxes with Borges et al. 2004
% parameterization (A = 0.27 cm/h)
%
% USAGE:-------------------------------------------------------------------  
% [Fd,F_wind,F_curr,k, k_wind, k_curr] = fas_B04(C,u10,uw,S,T,slp,gas,h,rh)
% [Fd,F_wind,F_curr,k, k_wind, k_curr] = fas_B04(0.25,7,2,35,20,1,'O2',1.5)
%   > Fd = -4.4860e-09
%   > Fc = 0
%   > Fp = 0
%   > Deq = 0
%   > k = 1.7365e-05
% 
% DESCRIPTION:-------------------------------------------------------------
% Calculate air-sea fluxes and steady-state supersaturation based on:
% Borges, A.V., Vanderborght, J., Schiettecatte, L. et al. Estuaries (2004) 
%   27: 593. doi:10.1007/BF02907647
%
% INPUTS:------------------------------------------------------------------
% 
% C:    gas concentration in mol m-3
% u10:  10 m wind speed (m/s)
% uw:   current speed (cm s-1)
% S:    Sea surface salinity
% T:    Sea surface temperature (deg C)
% slp:  sea level pressure (atm)
% gas:  formula for gas (He, Ne, Ar, Kr, Xe, N2, or O2), formatted as a
%       string, e.g. 'He'
% h:    water depth (m)
% rh:   relative humidity as a fraction of saturation (0.5 = 50% RH)
%       rh is an optional but recommended argument. If vs., it
%       will be automatically set to 1 (100% RH).
%
% OUTPUTS:-----------------------------------------------------------------
% Fd:       Diffusive air-sea flux                     (mol m-2 s-1)
% F_curr:   Flux from current derived turbulence       (mol m-2 s-1)
% F_wind:   Flux from wind derived turbulence          (mol m-2 s-1)
% 
% k:        Diffusive gas transfer velocity coef             (m s-1)
% k_curr:   Diffusive gas transfer velocity due to current   (m s-1)
% k_wind:   Diffusive gas transfer velocity due to wind      (m s-1)
% Note: Total air-sea flux is Fd = F_curr + F_wind
%
%
% REFERENCE:---------------------------------------------------------------
% Borges, A.V., Vanderborght, J., Schiettecatte, L. et al. Estuaries (2004) 
%   27: 593. doi:10.1007/BF02907647
%
% AUTHOR:---------------------------------------------------------------
% David Nicholson dnicholson@whoi.edu
% Woods Hole Oceanographic Institution
% Version: 12 April 2017
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2017 David Nicholson 
%
% Licensed under the Apache License, Version 2.0 (the "License"); you may 
% not use this file except in compliance with the License, which is
% available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fd,F_wind,F_curr,k, k_wind, k_curr] = fas_B04(C,u10,uw,S,T,slp,gas,h,rh)

% -------------------------------------------------------------------------
% Check for humidity, calculate dry pressure
% -------------------------------------------------------------------------

% if humidity is vs., set to 1 for all values
if nargin == 8
    rh = 1.0.*ones(size(C));
end;

% Equilibrium gas conc is referenced to 1 atm total air pressure, 
% including saturated water vapor (rh=1).
% Calculate ratio (observed dry air pressure)/(reference dry air pressure).
ph2oveq = vpress(S,T);
ph2ov = rh.*ph2oveq;
slpc = (slp-ph2ov)./(1-ph2oveq);

% -------------------------------------------------------------------------
% Calc gas exchange fluxes
% -------------------------------------------------------------------------

Geq = gasmoleq(S,T,gas);
[~,Sc] = gasmoldiff(S,T,gas);

cmh_2_ms = 100./(60*60);
k600_curr = cmh_2_ms.*1.719.*(uw./100).^-0.5.*(h.^-0.5);
k600_wind = cmh_2_ms.*(1.0 + 2.58.*u10);

n = 0.5;
k_curr = k600_curr.*(Sc./600).^-n;
k_wind = k600_wind.*(Sc./600).^-n;
k = k_curr + k_wind;

Fd = -k.*(C-Geq.*slpc);
F_curr = -k_curr.*(C-Geq.*slpc);
F_wind = -k_wind.*(C-Geq.*slpc);

