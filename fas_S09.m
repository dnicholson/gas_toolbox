% fas_S09: Function to calculate air-sea gas exchange flux using Stanley
% 2009 parameterization
%
% USAGE:-------------------------------------------------------------------
% [Fd, Fc, Fp, Deq, k] = fas_S09(C,u10,S,T,slp,gas,rh)
% [Fd, Fc, Fp, Deq, k] = fas_S09(0.01410,5,35,10,1,'Ar',0.9)
%   > Fd = -4.9960e-09
%   > Fc = 7.3493e-10
%   > Fp = 1.8653e-13
%   > Deq = 0.0027
%   > k = 1.9340e-05
%
% DESCRIPTION:-------------------------------------------------------------
% Calculate air-sea fluxes and steady-state supersaturation based on:
% Stanley, R.H., Jenkins, W.J., Lott, D.E., & Doney, S.C. (2009). Noble
% gas constraints on air-sea gas exchange and bubble fluxes. Journal of
% Geophysical Research: Oceans, 114(C11), doi: 10.1029/2009JC005396
%
% These estimates are valid over the range of wind speeds observed at
% Bermuda (0-13 m/s) and for open ocean, oligotrophic waters low in
% surfactants. Additionally, the estimates were determined using QuikSCAT
% winds, and if using another global wind product (e.g., NCEP reanalysis),
% a correction for biases between the wind products may be appropriate.
% Contact Rachel Stanley (rachel.stanley@wellesley.edu) with questions.
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
%
% INPUTS:------------------------------------------------------------------
% C:    gas concentration (mol/m^3)
% u10:  10 m wind speed (m/s)
% S:    Sea surface salinity (PSS)
% T:    Sea surface temperature (deg C)
% slp:  sea level pressure (atm)
% gas:  formula for gas (He, Ne, Ar, Kr, Xe, N2, or O2), formatted as a
%       string, e.g. 'He'
% rh:   relative humidity in the marine boundary layer as a fraction of
%       saturation (0.5 = 50% RH).
%       rh is an optional but recommended argument. If not provided, it
%       will be automatically set to 1 (100% RH).
%
% OUTPUTS:-----------------------------------------------------------------
% Fd:   Diffusive air-sea flux                        (mol m-2 s-1)
% Fc:   Flux from fully collapsing small bubbles      (mol m-2 s-1)
% Fp:   Flux from partially collapsing large bubbles  (mol m-2 s-1)
% Deq:  Equilibrium supersaturation                   (unitless (%sat/100))
% k:    Diffusive gas transfer velocity               (m s-1)
%
% Note: Total air-sea flux is Ft = Fd + Fp + Fc
%
% REFERENCE:---------------------------------------------------------------
%
% Stanley, R.H., Jenkins, W.J., Lott, D.E., & Doney, S.C. (2009). Noble
% gas constraints on air-sea gas exchange and bubble fluxes. Journal of
% Geophysical Research: Oceans, 114(C11), doi: 10.1029/2009JC005396
%
% Bubble penetration depth parameterization:
% Graham, A., D. K. Woolf, and A. J. Hall (2004), Aeration due to breaking
% waves. Part I: Bubble populations, J. Phys. Oceanogr., 34(5), 989?1007,
% doi:10.1175/1520-0485(2004)034<0989:ADTBWP>2.0.CO;2.
%
% AUTHOR:------------------------------------------------------------------
%
% Cara Manning (cmanning@whoi.edu) Woods Hole Oceanographic Institution
% Version: 12 April 2017
% Checked and approved by Rachel Stanley on September 20, 2015.
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2017 Cara Manning 
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fd, Fc, Fp, Deq, k] = fas_S09(C,u10,S,T,slp,gas,rh)
% -------------------------------------------------------------------------
% Conversion factors and constants
% -------------------------------------------------------------------------
atm2Pa = 1.01325e5; % Pascals per atm
R = 8.314;  % ideal gas constant in m3 Pa / (K mol)

% -------------------------------------------------------------------------
% Scaling factors for gas exchange coefficients
% -------------------------------------------------------------------------

Ac = 9.09E-11;
Ap = 2.29E-3;
gammaG = 0.97;
diffexp=2/3; betaexp=1;

% -------------------------------------------------------------------------
% Check for humidity
% -------------------------------------------------------------------------

% if humidity is not provided, set to 1 for all values
if nargin == 6
    rh =1.*ones(size(C));
end;

% -------------------------------------------------------------------------
% Calculate diffusive flux
% -------------------------------------------------------------------------

[D,Sc] = gasmoldiff(S,T,gas);
Geq = gasmoleq(S,T,gas);

k = gammaG*kgas(u10,Sc,'W92b'); % k_660 = 0.31 cm/hr

% slpc = (observed dry air pressure)/(reference dry air pressure)
% see Description section in header
ph2oveq = vpress(S,T);
ph2ov = rh.*ph2oveq;
slpc = (slp-ph2ov)./(1-ph2oveq);

% calculate diffusive flux with correction for local humidity
Fd = -k.*(C-Geq.*slpc);

% -------------------------------------------------------------------------
% Calculate complete trapping / air injection flux
% -------------------------------------------------------------------------

% air injection factor as a function of wind speed
% set to 0 below u10 = 2.27 m/s
wfact=(u10-2.27).^3; 
wfact(wfact<0) = 0;

% calculate dry atmospheric pressure in atm
patmdry=slp-ph2ov; % pressure of dry air in atm 

ai=gasmolfract(gas).*wfact.*patmdry.*atm2Pa./(R*(273.15+T));
Fc = Ac*ai; 


% -------------------------------------------------------------------------
% Calculate partial trapping / exchange flux
% -------------------------------------------------------------------------

% calculate bubble penetration depth, Zbub, then calculate hydrostatic
% pressure in atm
Zbub = 0.15*u10 - 0.55;
Zbub(Zbub<0)=0; 
phydro=(gsw_sigma0(S,T)+1000).*9.81.*Zbub./atm2Pa; 

% multiply by scaling factor Ap by beta raised to power betaexp and 
% diffusivity raised to power diffexp
apflux=ai.*Ap.*D.^diffexp.*(gasBunsen(S,T,gas).^betaexp); 
Fp=apflux.*(phydro./patmdry-C./Geq+1); 


% -------------------------------------------------------------------------
% Calculate steady-state supersaturation
% -------------------------------------------------------------------------
Deq = ((Fc+Fp)./k)./Geq;

end
