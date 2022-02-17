function [conc_CH4] = CH4sol(S,T,varargin)

% CH4sol   Solubility of CH4 in sea water
%=========================================================================
% CH4sol Version 1.1 3/5/2009
%          
% USAGE:  [conc_CH4]=CH4sol(S,T,varargin)
%
% DESCRIPTION:
%    Solubility (saturation) of methane (CH4) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%   (OPTIONAL:)
%   varargin = fCH4, optional input allowing user to specify fugacity of CH4 [atm]
% 
% DEFAULT, if fCH4 is not specified:
%   fCH4 = 1940e-9 atm
%
% OUTPUT:
%   concCH4 = solubility of CH4  [umol/kg] 
%
% Atm. concentration available at: NOAA / ESRL interactive data viewer
% (https://www.esrl.noaa.gov/gmd/dv/iadv/)
% 
% AUTHOR:  Cara Manning and P. Tortell modified from Roberta Hamme
%
% REFERENCE:
%    Denis Wiesenburg and Norman Guinasso, 1979.
%    "Equilibrium solubilities of methane, carbon monoxide and hydrogen in water and seawater"
%    Journal of chemical and engineering data, 24(4), pp. 356-360.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER: general purpose
% CALLEE: none

%----------------------
% Check input parameters
%----------------------
if nargin <2
   error('CH4sol.m: Must pass 2 or 3 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) || (ns~=nt)) && (ms+ns>2) && (mt+nt>2)
   error('CH4sol: S & T must have same dimensions or be singular')
end %if


% if CH4 fugacity is not provided, use default value
if nargin > 2
    fCH4 = varargin{1};
else
    fCH4 = 1940e-9;
end

%------
% BEGIN
%------

% convert T to absolute temperature
T=T+273.15;

% constants from Table VI of  (1979)  units in nmol/kg

A1 = -417.5053;
A2 = 599.8626;
A3 = 380.3636;
A4 = -62.0764;
B1 = -0.064236;
B2 = 0.034980;
B3 = -0.0052732;

%Fg=0.675e-006;
% Eqn (7) of Wiesenburg and Guinasso 1979
for i=1:length(T)
    if numel(fCH4)>1
    conc_CH4(i) = log(fCH4(i)) + A1 + A2*(100/T(i)) + A3*log(T(i)/100) + A4*(T(i)/100) + S(i)*(B1 + B2*(T(i)/100) + B3*(T(i)/100)^2);
    conc_CH4(i)=exp(conc_CH4(i));
    else
    conc_CH4(i) = log(fCH4) + A1 + A2*(100/T(i)) + A3*log(T(i)/100) + A4*(T(i)/100) + S(i)*(B1 + B2*(T(i)/100) + B3*(T(i)/100)^2);
    conc_CH4(i)=exp(conc_CH4(i));
    end
end

% convert from nmol/kg to umol/kg
conc_CH4 = conc_CH4./1000;

return