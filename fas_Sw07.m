% [Fd, Fc, Fp, Deq, k] = fas_Sw07(C,u10,S,T,slp,gas,rh)
% Function to calculate air-sea gas exchange flux using Sweeney 07
% parameterization (A = 0.27 cm/hr)
%
% USAGE:-------------------------------------------------------------------  
% [Fd, Fc, Fp, Deq, k] = fas_Sw07(C,u10,S,T,slp,gas,rh)
% [Fd, Fc, Fp, Deq, k] = fas_Sw07(0.01410,5,35,10,1,'Ar',0.9)
%   > Fd = -4.4860e-09
%   > Fc = 0
%   > Fp = 0
%   > Deq = 0
%   > k = 1.7365e-05
% 
% DESCRIPTION:-------------------------------------------------------------
% Calculate air-sea fluxes and steady-state supersaturation based on:
% Sweeney, C., Gloor, E., Jacobson, A. R., Key, R. M., McKinley, G.,
% Sarmiento, J. L., & Wanninkhof, R. (2007). Constraining global air?sea
% gas exchange for CO2 with recent bomb 14C measurements. Global
% Biogeochemical Cycles, 21(2).
%
% INPUTS:------------------------------------------------------------------
% 
% C:    gas concentration in mol m-3
% u10:  10 m wind speed (m/s)
% S:    Sea surface salinity
% T:    Sea surface temperature (deg C)
% slp:  sea level pressure (atm)
% gas:  formula for gas (He, Ne, Ar, Kr, Xe, N2, or O2), formatted as a
%       string, e.g. 'He'
% rh:   relative humidity as a fraction of saturation (0.5 = 50% RH)
%       rh is an optional but recommended argument. If not provided, it
%       will be automatically set to 0.8.
%
% OUTPUTS:-----------------------------------------------------------------
% Fd:   Diffusive air-sea flux                        (mol m-2 s-1)
% Fc:   Flux from fully collapsing small bubbles      (mol m-2 s-1)
% Fp:   Flux from partially collapsing large bubbles  (mol m-2 s-1)
% Deq:  Equilibrium supersaturation                   (unitless (%sat/100))
% k:    Diffusive gas transfer velocity               (m s-1)
%
% Note: Total air-sea flux is Ft = Fd + Fc + Fp
%
% Note: Fp, Fc, and Deq will always be 0 for fas_Sw07. They are
% included as outputs for consistency with the other fas functions.
%
% REFERENCE:---------------------------------------------------------------
% Sweeney, C., Gloor, E., Jacobson, A. R., Key, R. M., McKinley, G.,
% Sarmiento, J. L., & Wanninkhof, R. (2007). Constraining global air?sea
% gas exchange for CO2 with recent bomb 14C measurements. Global
% Biogeochemical Cycles, 21(2).
%
% AUTHOR:---------------------------------------------------------------
% Cara Manning cmanning@whoi.edu
% Woods Hole Oceanographic Institution
% Version: August 2015
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2015 Cara Manning 
%
% Licensed under the Apache License, Version 2.0 (the "License"); you may 
% not use this file except in compliance with the License, which is
% available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fd, Fc, Fp, Deq, k] = fas_Sw07(C,u10,S,T,slp,gas,rh)

% -------------------------------------------------------------------------
% Check for humidity, calculate dry pressure
% -------------------------------------------------------------------------

% if humidity is not provided, set to 0.8 for all values
if nargin == 6
    rh = 0.8.*ones(size(C));
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
k = kgas(u10,Sc,'Sw07');
Fd = -k.*(C-Geq.*slpc);

% Set Fc and Fp to 0. They are included to be consistent with the other
% fas_ functions.
Fc = zeros(size(Fd));
Fp = zeros(size(Fd));

% Calculate steady state equilibrium supersaturation. This will also be 0.
Deq = ((Fc+Fp)./k)./Geq;