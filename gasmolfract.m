% =========================================================================
% gasmolfract
% USAGE:
% Xg = gasmolfract(gas)   OR
% Xg = gasmolfract(gas,pCO2)
% -------------------------------------------------------------------------
% mole fraction in dry atmosphere of a well mixed atmospheric gas
% 
% Noble gases (excluding Ar) from:
% Glueckauf, E. (1951). The composition of atmospheric air. Compendium of
% Meteorology. 3-10.
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
% Helium:       'He'
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
% Xg        volumetric mixing of the gas in dry atm (mixing ratio)
%           equivalent to mol/mol if you assume gases are ideal
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% X = gasmolfract('Ar')
% X = 0.009332
%
% written by Roo Nicholson 08/03/08
% Also see: gasmoleq.m, gasmolsol.m
% =========================================================================

function [Xg] = gasmolfract(gas,varargin)

if nargin == 2
    pCO2 = 1e-6.*varargin{1};
elseif nargin == 1
    pCO2 = 1e-6.*400;
end

% these values are appropriate for pCO2 = 400 uatm
if strcmpi(gas, 'He')
    Xg = 5.24e-6;
elseif strcmpi(gas, 'Ne')
    Xg = 0.00001818;
elseif strcmpi(gas, 'Ar')
    Xg = 0.009332;
elseif strcmpi(gas, 'Kr')
    Xg = 0.00000114;
elseif strcmpi(gas, 'Xe')
    Xg = 8.7e-8;
elseif strcmpi(gas, 'N2')
    Xg = 0.780848;
elseif strcmpi(gas, 'Ar36')
    Xg = 0.009332.*0.003651267;
elseif strcmpi(gas, 'O2')
    Xg = 0.209790-pCO2;
else
    error('Gas name must be He, Ne, Ar, Kr, Xe, N2, O2 or Ar36');
end
