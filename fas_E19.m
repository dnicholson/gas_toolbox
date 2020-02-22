% fas_E19: Function to calculate air-sea fluxes with Liang 2013
% parameterization
%
% USAGE:-------------------------------------------------------------------
% [Fd, Fc, Fp, Deq, k] = fas_E19(C,u10,S,T,slp,gas,rh)
% [Fd, Fc, Fp, Deq, k] = fas_E19(0.01410,5,35,10,1,'Ar')
%   >Fd = -5.6030e-09
%   >Fc = 5.0339e-11
%   >Fp = -2.4485e-10
%   >Deq = 0.0012
%   >k = 2.0377e-05
%
%
% DESCRIPTION:-------------------------------------------------------------
%
% Calculate air-sea fluxes and steady-state supersat based on Emerson et
% al. (2019) modification to the Liang et al. 2013 parameterization which
% multiplies bubble fluxes by 0.37
%
%
% INPUTS:------------------------------------------------------------------
% C:    gas concentration (mol m-3)
% u10:  10 m wind speed (m/s)
% SP:   Sea surface salinity (PSS)
% pt:   Sea surface temperature (deg C)
% pslp: sea level pressure (atm)
% gas:  formula for gas (He, Ne, Ar, Kr, Xe, N2, or O2), formatted as a
%       string, e.g. 'He'
% rh:   relative humidity as a fraction of saturation (0.5 = 50% RH)
%       rh is an optional but recommended argument. If not provided, it
%       will be automatically set to 1 (100% RH).
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
% OUTPUTS:-----------------------------------------------------------------
%
% Fd:   Surface gas flux                              (mol m-2 s-1)
% Fc:   Flux from fully collapsing small bubbles      (mol m-2 s-1)
% Fp:   Flux from partially collapsing large bubbles  (mol m-2 s-1)
% Deq:  Equilibrium supersaturation                   (unitless (%sat/100))
% k:    Diffusive gas transfer velocity               (m s-1)
%
% Note: Total air-sea flux is Ft = Fd + Fc + Fp
%
% REFERENCE:---------------------------------------------------------------
% Emerson, S., Yang, B., White, M., & Cronin, M. (2019). Air?sea gas 
%   transfer: Determining bubble fluxes with in situ N2 observations. 
%   Journal of Geophysical Research: Oceans, 124, 2716?2727. 
%   https://doi.org/10.1029/2018JC014786
%
% Liang, J.-H., C. Deutsch, J. C. McWilliams, B. Baschek, P. P. Sullivan, 
%   and D. Chiba (2013), Parameterizing bubble-mediated air-sea gas 
%   exchange and its effect on ocean ventilation, Global Biogeochem. Cycles, 
%   27, 894?905, doi:10.1002/gbc.20080.
%
% AUTHOR:---------------------------------------------------------------
% Written by David Nicholson dnicholson@whoi.edu
%
% Woods Hole Oceanographic Institution
% Version: 21 Feb 2020
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2017 David Nicholson and Cara Manning 
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fd, Fc, Fp, Deq, Ks] = fas_E19(C,u10,SP,pt,pslp,gas,rh)

% Bubble scaling factor from Emerson et al. 2019
bfact = 0.37;

% -------------------------------------------------------------------------
% Conversion factors
% -------------------------------------------------------------------------
m2cm = 100; % cm in a meter
h2s = 3600; % sec in hour
atm2Pa = 1.01325e5; % Pascals per atm

% -------------------------------------------------------------------------
% Calculate water vapor pressure and adjust sea level pressure
% -------------------------------------------------------------------------

% if humidity is not provided, set to 1 for all values
if nargin == 6
    rh =ones(size(C));
end

ph2oveq = vpress(SP,pt);
ph2ov = rh.*ph2oveq;

% slpc = (observed dry air pressure)/(reference dry air pressure)
% see Description section in header of fas_N11.m
pslpc = (pslp - ph2ov)./(1 - ph2oveq);

% -------------------------------------------------------------------------
% Parameters for COARE 3.0 calculation
% -------------------------------------------------------------------------

% Calculate potential density at surface
SA = SP.*35.16504./35;
CT = gsw_CT_from_pt(SA,pt);
rhow = gsw_sigma0(SA,CT)+1000;
rhoa = 1.225;

lam = 13.3;
A = 1.3;
phi = 1;
tkt = 0.01;
hw=lam./A./phi;
ha=lam;

% air-side schmidt number
ScA = 0.9;

R = 8.314;  % units: m3 Pa K-1 mol-1

% -------------------------------------------------------------------------
% Calculate gas physical properties
% -------------------------------------------------------------------------
xG = gasmolfract(gas);
Geq = gasmoleq(SP,pt,gas);
alc = (Geq/atm2Pa).*R.*(pt+273.15);

Gsat = C./Geq;
[~, ScW] = gasmoldiff(SP,pt,gas);

% -------------------------------------------------------------------------
% Calculate COARE 3.0 and gas transfer velocities
% -------------------------------------------------------------------------
% ustar
cd10 = cdlp81(u10);
ustar = u10.*sqrt(cd10);

% water-side ustar
ustarw = ustar./sqrt(rhow./rhoa);

% water-side resistance to transfer
rwt = sqrt(rhow./rhoa).*(hw.*sqrt(ScW)+(log(.5./tkt)/.4));

% air-side resistance to transfer
rat = ha.*sqrt(ScA)+1./sqrt(cd10)-5+.5*log(ScA)/.4;

% diffusive gas transfer coefficient (L13 eqn 9)
Ks = ustar./(rwt+rat.*alc);

% bubble transfer velocity (L13 eqn 14)
Kb = bfact.*1.98e6.*ustarw.^2.76.*(ScW./660).^(-2/3)./(m2cm.*h2s);

% overpressure dependence on windspeed (L13 eqn 16)
dP = 1.5244.*ustarw.^1.06;


% -------------------------------------------------------------------------
% Calculate air-sea fluxes
% -------------------------------------------------------------------------

Fd = Ks.*Geq.*(pslpc-Gsat); % Fs in L13 eqn 3
Fp = Kb.*Geq.*((1+dP).*pslpc-Gsat); % Fp in L13 eqn 3
Fc = bfact.*xG.*5.56.*ustarw.^3.86; % L13 eqn 15

% -------------------------------------------------------------------------
% Calculate steady-state supersaturation 
% -------------------------------------------------------------------------
Deq = (Kb.*Geq.*dP.*pslpc+Fc)./((Kb+Ks).*Geq.*pslpc); % L13 eqn 5

end

function [ cd ] = cdlp81( u10)
% Calculates drag coefficient from u10, wind speed at 10 m height

cd = (4.9e-4 + 6.5e-5 * u10);
cd(u10 <= 11) = 0.0012;
cd(u10 >= 20) = 0.0018;

end