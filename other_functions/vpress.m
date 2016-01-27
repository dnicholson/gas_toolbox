function [vapor_press_atm] = vpress(S,T)

% vpress   Vapor pressure of sea water
%=========================================================================
% vpress Version 2.0 : 27 October 2012
%          Author: Roberta C. Hamme (University of Victoria)
%
% USAGE:  vapor_press = vpress(S,T)
%
% DESCRIPTION:
%    Vapor pressure of sea water
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS-78]
%   T = temperature [degree C]
%
% OUTPUT:
%   vapor_press = vapor pressure of seawater  [atm] 
% 
% AUTHOR:  Roberta Hamme (rhamme@uvic.ca)
%
% REFERENCE:
%   Guide to Best Practices for Ocean CO2 Measurements
%   Dickson, A.G., C.L. Sabine, J.R. Christian (Eds.) 2007
%   PICES Special Publication 3, 191pp.
%   Chapter 5: Physical and thermodynamic data
%       Based on: Wagner, W., A. Pruss (2002) The IAPWS formulation 1995 
%       for the thermodynamic properties of ordinary water substance for 
%       general and scientific use, J. Phs. Chem. Ref. Data, 31, 387-535.
%       AND Millero, F.J. (1974) Seawater as a multicomponent electrolyte 
%       solution, pp.3-80.  In: The Sea, Vol. 5, E.D. Goldberg Ed.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER: general purpose
% CALLEE: none

%----------------------
% Check input parameters
%----------------------
if nargin ~=2
   error('vpress.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
% Check that T&S have the same shape or are singular
if ((ms~=mt) || (ns~=nt)) && (ms+ns>2) && (mt+nt>2)
   error('vpress: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

%Calculate temperature in Kelvin and modified temperature for Chebyshev polynomial
temp_K = T+273.15;
temp_mod = 1-temp_K./647.096;

%Calculate value of Wagner polynomial
Wagner = -7.85951783*temp_mod +1.84408259*temp_mod.^1.5 -11.7866497*temp_mod.^3 +22.6807411*temp_mod.^3.5 -15.9618719*temp_mod.^4 +1.80122502*temp_mod.^7.5;

%Vapor pressure of pure water in kiloPascals and mm of Hg
vapor_0sal_kPa = exp(Wagner * 647.096 ./ temp_K) .* 22.064 * 1000;

%Correct vapor pressure for salinity
molality = 31.998 * S ./(1e3-1.005*S);
osmotic_coef = 0.90799 -0.08992*(0.5*molality) +0.18458*(0.5*molality).^2 -0.07395*(0.5*molality).^3 -0.00221*(0.5*molality).^4;
vapor_press_kPa = vapor_0sal_kPa .* exp(-0.018 * osmotic_coef .* molality);

%Convert to atm
vapor_press_atm = vapor_press_kPa/101.32501;