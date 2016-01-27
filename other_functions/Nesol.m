function [conc_Ne] = Nesol(S,T)

% Nesol   Solubility of Ne in sea water
%=========================================================================
% Nesol Version 1.1 4/4/2005
%          Author: Roberta C. Hamme (Scripps Inst of Oceanography)
%
% USAGE:  concNe = Nesol(S,T)
%
% DESCRIPTION:
%    Solubility (saturation) of neon (Ne) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%
% OUTPUT:
%   concNe = solubility of Ne  [umol/kg] 
% 
% AUTHOR:  Roberta Hamme (rhamme@ucsd.edu)
%
% REFERENCE:
%    Roberta Hamme and Steve Emerson, 2004.
%    "The solubility of neon, nitrogen and argon in distilled water and seawater."
%    Deep-Sea Research I, 51(11), p. 1517-1528.
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
   error('Nesol.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('Nesol: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

% convert T to scaled temperature
temp_S = log((298.15 - T)./(273.15 + T));

% constants from Table 4 of Hamme and Emerson 2004
A0_neon = 2.18156;
A1_neon = 1.29108;
A2_neon = 2.12504;
B0_neon = -5.94737e-3;
B1_neon = -5.13896e-3;

% Eqn (1) of Hamme and Emerson 2004
conc_Ne = exp(A0_neon + A1_neon*temp_S + A2_neon*temp_S.^2 + S.*(B0_neon + B1_neon*temp_S));

% Convert from nmol/kg to umol/kg
conc_Ne = conc_Ne/1000;

return