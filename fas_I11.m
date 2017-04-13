% fas_N11: Function to calculate air-sea gas exchange flux using Nicholson
% 2011 parameterization
%
% USAGE:-------------------------------------------------------------------
%  
% [Fd, Fc, Fp, Deq, k] = fas_I11(C,u10,S,T,slp,gas,rh)
% [Fd, Fc, Fp, Deq, k] = fas_I11(0.01410,5,35,10,1,'Ar',1)
%
% > Fd = -5.9209e-09
% > Fc = 4.5048e-09
% > Fp = 0
% > Deq = 0.0151
% > k = 2.1533e-05
%
% DESCRIPTION:-------------------------------------------------------------
%
% Calculate air-sea fluxes and steady-state supersaturation based on:
% Ito, T., R. C. Hamme, and S. Emerson (2011), Temporal and spatial 
%   variability of noble gas tracers in the North Pacific, J. Geophys. 
%   Res., 116, C08039, doi:10.1029/2010JC006828.
%
% Fc = Ainj * slpc * Xg * u3
% Fp = 0
%
% where u3 = ((u-2.27)./(10-2.27))^3 (and zero for  u < 2.27)
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
% 
% C:    gas concentration in mol m-3
% u10:  10 m wind speed (m/s)
% S:    Sea surface salinity (PSS)
% T:    Sea surface temperature (deg C)
% slp:  sea level pressure (atm)
%
% gas:  formula for gas, formatted as a string, e.g. 'He'
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
% varargin: optional but recommended arguments
%       rhum: relative humidity as a fraction of saturation (0.5 = 50% RH).
%             If not provided, it will be automatically set to 1 (100% RH).
%
% OUTPUTS:-----------------------------------------------------------------
% Fd:   Surface air-sea diffusive flux                [mol m-2 s-1]
% Fc:   Injection bubble flux (complete trapping)     [mol m-2 s-1]
% Fp:   Exchange bubble flux (partial trapping)       [mol m-2 s-1]
% Deq:  Steady-state supersaturation                  [unitless (%sat/100)]
% k:    Diffusive gas transfer velocity               (m s-1)
%
% Note: Total air-sea flux is Ft = Fd + Fc + Fp
%
% REFERENCE:---------------------------------------------------------------
% Ito, T., R. C. Hamme, and S. Emerson (2011), Temporal and spatial 
%   variability of noble gas tracers in the North Pacific, J. Geophys. 
%   Res., 116, C08039, doi:10.1029/2010JC006828.
%
% AUTHORS:-----------------------------------------------------------------
% David Nicholson dnicholson@whoi.edu
% Woods Hole Oceanographic Institution
% Version 12 April 2017
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2017 David Nicholson  
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fd, Fc, Fp, Deq, k] = fas_I11(C,u10,S,T,slp,gas,varargin)

% see Appendix A in Ito et al. 2011
Ainj = 9.1e-9; % m s-1
alph = 9.3e-5; % m s-1
Aex = 0;
R = 8.314;
Patm = 1.01325e5; % 1 atm in Pa

% if humidity is not provided, set to 1 for all values
if nargin > 6
    rhum = varargin{1};
else
    rhum = 1;
end

% slpc = (observed dry air pressure)/(reference dry air pressure)
% see Description section in header
ph2oveq = vpress(S,T);
ph2ov = rhum.*ph2oveq;
slpc = (slp-ph2ov)./(1-ph2oveq);
slpd = slp-ph2ov;

[D,Sc] = gasmoldiff(S,T,gas);
Geq = gasmoleq(S,T,gas);

% calculate wind speed term for bubble flux
u3 = ((u10-2.27)./(10-2.27)).^3;
u3(u3 < 0) = 0;

k = alph.*(u10./10).^2.*(Sc./660).^-0.5;
Fd = -k.*(C-slpc.*Geq);
Fc = (Ainj./(R.*(T+273.15))).*slpd.*Patm.*gasmolfract(gas).*u3;
Fp = Aex.*slpc.*Geq.*D.^0.5.*u3;
Deq = ((Fp+Fc)./k)./Geq;

