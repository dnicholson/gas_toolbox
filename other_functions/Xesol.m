function [conc_Xe] = Xesol(S,T)

% Xesol   Solubility of Xe in sea water
%=========================================================================
% Xesol Version 1.2 1/11/2007
%          Author: Roberta C. Hamme (University of Victoria)
%
% USAGE:  concXe = Xesol(S,T)
%
% DESCRIPTION:
%    Solubility (saturation) of xenon (Xe) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%
% OUTPUT:
%   concXe = solubility of Xe  [umol/kg] 
% 
% AUTHOR:  Roberta Hamme (rhamme@ucsd.edu)
%
% REFERENCE:
%    R. Hamme fit to data of
%    D. Wood and R. Caputi (1966) "Solubilities of Kr and Xe in fresh and sea water"
%    U.S. Naval Radiological Defense Laboratory, Technical Report USNRDL-TR-988, 
%    San Francisco, CA, pp. 14.
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
   error('Xesol.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('Xesol: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

% convert T to scaled temperature
temp_S = log((298.15 - T)./(273.15 + T));

% constants from fit procedure of Hamme and Emerson 2004 to Wood and Caputi data
A0_xenon = -7.48588;
A1_xenon = 5.08763;
A2_xenon = 4.22078;
B0_xenon = -8.17791e-3;
B1_xenon = -1.20172e-2;

% Eqn (1) of Hamme and Emerson 2004
conc_Xe = exp(A0_xenon + A1_xenon*temp_S + A2_xenon*temp_S.^2 + S.*(B0_xenon + B1_xenon*temp_S));

return