% =========================================================================
% [pGdry] = gasvolfract(gas,varargin)
% -------------------------------------------------------------------------
% volume ratio in dry atmosphere of a well mixed atmospheric gas
% 
% Noble gases (excluding Ar) from:
% Glueckauf, E. (1951) The composition of atmospheric air.
% T.F. Malone (Ed.), Compendium of Meteorology, Amer. Meteorological Soc,
% Boston (1951), pp. 3?11
%
% N2, Ar, O2 are from Table 1 in:
% Picard, A., Davis, R. S., Gläser, M., & Fujii, K. (2008). Revised formula 
%   for the density of moist air (CIPM-2007). Metrologia, 45(2), 149.
%   doi:10.1088/0026-1394/45/2/004
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% gas:      name of the gas (see below)
%
% Neon:         'Ne'
% Argon:        'Ar'
% Krypton:      'Kr'
% Xenon:        'Xe'
% Nitrogen:     'N2'
% Oxygen:       'O2'
%
% optional input:
% pCO2      pCO2 in uatm (or ppmv)  This is used to calculate dilution due
%           to addded CO2 - assumes xCO2 + xO2 are a constant
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% pGdry     volumetric mixing of the gas in dry atm (mixing ratio)
%           this uses non-ideal molar volumes from gasmolvol.m
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% X = gasmolefract('Ar')
% X = 0.00934
%
% written by Roo Nicholson 16 SEP 2015 dnicholson@whoi.edu
% Also see: gasmolvol.m, gasmoleq.m, gasmolsol.m
% =========================================================================


function [pGdry] = gasvolfract(gas,varargin)


pCO2ref = 1e-6.*400;
if nargin == 2
    pCO2 = 1e-6.*varargin{1};
elseif nargin == 1
    pCO2 = pCO2ref;
end

% below is code used to calculate 'vSum' for pCO2 = 400 uatm
% commmented out to avoid recomputing for every call

% vN2 = gasmolfract('N2').*gasmolvol('N2');
% vO2 = gasmolfract('O2').*gasmolvol('O2');
% vHe = gasmolfract('He').*gasmolvol('He');
% vNe = gasmolfract('Ne').*gasmolvol('Ne');
% vAr = gasmolfract('Ar').*gasmolvol('Ar');
% vKr = gasmolfract('Kr').*gasmolvol('Kr');
% vXe = gasmolfract('Xe').*gasmolvol('Xe');
% vCO2 = pCO2.*gasmolvol('CO2');
% 
% xsum = gasmolefract('N2')+gasmolefract('O2')+gasmolefract('Ar')+...
%     gasmolefract('He')+gasmolefract('Ne')+gasmolefract('Kr')+...
%     gasmolefract('Xe');
% 
% vSum = (vN2+vO2+vHe+vNe+vAr+vKr+vXe+vCO2)./xsum;

% average molar volume of air
vSum = 22.410316382164048;
if nargin == 2
    % this is a tiny correction 
    vSum = vSum + (pCO2-pCO2ref).*(gasmolvol('CO2')-gasmolvol('O2'));
end
switch gas
    case {'He','Ne','Ar','Kr','Xe','N2','O2'}
        pGdry = gasmolfract(gas).*gasmolvol(gas)./vSum;
    otherwise
    error('Gas name must be He, Ne, Ar, Kr, Xe, N2, O2');
end
