function [conc_n2o, Kh_n2o] = N2Osol(S,T,fN2O)

% N2O sol   Solubility of nitrous oxide in sea water
%=========================================================================
% 
% USAGE:  [conc_n2o, Kh_n2o]=n2osol(S,T,fn2o)
%   fn2o = optional
%
% DESCRIPTION:
%    Solubility (saturation) of nitrous oxide (N2O) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%   (OPTIONAL:)
%   fn2o = allows user to specify the N2O fugacity [atm]
% 
% DEFAULT, if fN2O is not specified: 
%   fN2O = 329e-9 atm
%
% OUTPUT:
%   conc_n2o = solubility of N2O  [umol/kg] 
%   Kh_n2o = Henry's law constant for N2O [mol kg^-1 atm^-1]
% 
% AUTHOR:  P. Tortell modified from Roberta Hamme
% Edited by: C. Manning to convert units to umol/kg
% 
% Atm. concentration available at: NOAA / ESRL interactive data viewer
% (https://www.esrl.noaa.gov/gmd/dv/iadv/)
% 
% REFERENCE:
%    Marine Chemistry,8(1980), 347-359 
%    "NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER "
%    Ray Weiss and B.A. Price
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER: general purpose
% CALLEE: none

%----------------------
% Check input parameters
%----------------------
% if nargin ~=2
if nargin <2
   error('N2Osol.m: Must pass at least 2 parameters')
end %if

% Check S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
  
% Check that T&S have the same shape or are singular
if ((ms~=mt) || (ns~=nt)) && (ms+ns>2) && (mt+nt>2)
   error('N2Osol: S & T must have same dimensions or be singular')
end 

% If fN2O is not specified use default value
if nargin==2
    fN2O = 329e-9; % 329 parts per billion
end

%------------------
% BEGIN CALCULATION
%------------------

% convert T to absolute temperature
T=T+273.15;

% constants from Table II of Weiss and Price units in mol kg^-1 atm^-1

A1 = -64.8539;
A2 = 100.2520;
A3 = 25.2049;
B1 = -0.062544;
B2 = 0.035337;
B3 = -0.0054699;

% Eqn (7) of Wiesenburg and Guinasso 1979
Kh_n2o = nan.*T;
for i=1:length(T)
Kh_n2o(i) =  A1 + A2*(100/T(i)) + A3*log(T(i)/100) + S(i)*(B1 + B2*(T(i)/100) + B3*(T(i)/100)^2);
Kh_n2o(i)=exp(Kh_n2o(i));
end

conc_n2o=Kh_n2o.*fN2O.*1e6;

return