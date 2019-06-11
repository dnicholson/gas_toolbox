% [D, Sc] = gasmoldiff(S,T,gas)
% Diffusion coeff and Schmidt number for gases in fresh/sea water
%=========================================================================
% Modified from gas_diffusion Version 2.0 16 July 2013
%          Author: Roberta C. Hamme (University of Victoria)
%
% USAGE:  [DAr ScAr] = gasmoldiff(35,20,'Ar')
% 
% > DAr = 2.2609e-09
% > ScAr = 463.1527
%
%
% DESCRIPTION:
%    Diffusion coefficients of various gases in fresh/sea water
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS-78]
%   T = temperature [degree C]
%   gas = 'He','He','Ne','Ar','Kr','Xe','N2','O2','CH4','N2O',or 'CO2'
% 
% OUTPUT:
%   D = diffusion coefficient   [m^2/s] 
%   Sc = Schmidt number
% 
% AUTHOR:  David Nicholson : dnicholson@whoi.edu
%
% REFERENCE:
%    He, Ne, Kr, Xe, CH4 freshwater values from Jahne et al., 1987.
%       "Measurement of Diffusion Coeffients of Sparingly Soluble Gases in Water"
%       J. Geophys. Res., 92(C10), 10767-10776.
%    Ar freshwaters values are extrapolated from Jahne et al. 1987
%       He, Ne, Kr, Xe values at each temperature were fitted to D vs. mass^-0.5 
%       relationship to predict Ar at those temperatures, then Ar was fit to a 
%       ln(D_Ar) vs. 1/T(K) relationship to obtain Eyring equation coefficients
%    O2 and N2 freshwater values from  Ferrell and Himmelblau, 1967.
%       "Diffusion coefficients of nitrogen and oxygen in water"
%       J. Chem. Eng. Data, 12(1), 111-115, doi: 10.1021/je60032a036.
%    Correction for salinity is based on Jahne's observed average 4.9% decrease in 
%       diffusivity for H2 and He in 35.5 ppt NaCl solution
%
%    for Ne, the Jahne values compare well with and fall between those of
%       Wise and Houghton 1968 and Holz et al. 1994
%    for Ar, the extrapolated Jahne values compare well with Wise and Houghton 1968,
%       O'Brien and Hyslop 1977, and a numerical simulation by Bourg et al. 2008
%       but are higher than other reported values
%    for Kr, the Jahne values compare well with Wise and Houghton 1968,
%       and a numerical simulation by Bourg et al. 2008
%    for Xe, the Jahne values compare well with Pollack 1981, and a numerical 
%       simulation by Bourg et al. 2008, but fall significantly above Wise and Houghton 1968
%       and below Weingartner et al. 1992
%    for O2, there is general agreement among measurements. The Ferrel and Himmelblau values 
%       agree reasonably well with Baird and Davidson 1962, Wise and Houghton 1966,
%       Duda and Vrentas 1968, O'Brien and Hyslop 1977, and the Wilke and Change (1955) theory 
%       as tabulated by Wanninkhof 1992, but lie below Krieger et al 1967
%    for N2, there is less agreement. The Ferrel and Himmelblau values 
%       agree reasonably well with Baird and Davidson 1962, O'Brien and Hyslop 1977, 
%       and the Wilke and Change (1955) theory as tabulated by Wanninkhof 1992, 
%       but lie significantly below the values of Wise and Houghton 1966 and Krieger et al 1967
%    for He, I did not investigate comparisons of data, but chose Jahne 
%       since their work for other gases appears to be the best
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

function [D, Sc] = gasmoldiff(S,T,gas)

R = 8.314510;

if strcmpi(gas, 'He')
    [AEa] = [0.8180e-6 11700];
elseif strcmpi(gas, 'Ne')
    [AEa] = [1.6080e-6 14840];
elseif strcmpi(gas, 'O2')
    [AEa] = [4.286e-6 18700];
elseif strcmpi(gas, 'Ar')
    [AEa] = [2.227e-6 16680];
elseif strcmpi(gas, 'Ar36')
    [AEa] = [2.227e-6 16680];
elseif strcmpi(gas, 'Kr')
    [AEa] = [6.3930e-6 20200];
elseif strcmpi(gas, 'Xe')
    [AEa] = [9.0070e-6 21610];
elseif strcmpi(gas, 'N2')
    [AEa] = [3.4120e-6 18500];
elseif strcmpi(gas, 'CH4')
    [AEa] = [3.0470e-6 18360];
elseif strcmpi(gas, 'H2')
    [AEa] = [3.3380e-6 16060];
elseif strcmpi(gas, 'CO2')
    % CO2 schmidt# formula from Wanninkhof (1992)
    Sc = 2073.1 - 125.62.*T + 3.6276.*T.^2 - 0.043219.*T.^3;
    D = sw_visc(S,T,0)./Sc;
    return
    
elseif strcmpi(gas,'N2O')
    % Wanninkhof (2014) 4th order polynomial
    % fit using Wilke and Chang (1955) as adapted by Hayduk and Laudie (1974)
    a = [2356.2; 166.38; 6.3952; 0.13422; 0.0011506];
    Sc = a(1) - a(2) .* T + a(3) .* T.^2 - a(4) .* T.^3 + a(5) .* T.^4;
    D = sw_visc(S,T,0)./Sc;
    return
else
    error('Gas name must be He, Ne, O2, N2, Ar, Kr, Xe, CO2, CH4, N2O, or H2');
end

%freshwater diffusivity
D0 = AEa(1).*exp(-AEa(2)./(R.*(T+273.16)));
D = D0 .* (1 - 0.049 * S / 35.5);

Sc = sw_visc(S,T,0)./D;
if strcmpi(gas, 'Ar36')
    alphaD = 0.995;
    D = D./alphaD^2;
    Sc = Sc.*alphaD^2;
end


