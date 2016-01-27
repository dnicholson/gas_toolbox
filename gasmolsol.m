% =========================================================================
% GASMOLSOL.M - calculates Henry's Law solubility (for a pure gas)
% in mol m-3 atm-1 
%
% This is a wrapper function. See individual solubility functions for more
% details.
%
% [sol] = gasmolsol(SP,pt,gas)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% SP        Practical Salinity
% pt        Potential temperature [degC]
% gas       gas string: He, Ne, Ar, Kr, Xe, O2 or N2
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sol       Henry's Law solubility in mol m-3 atm-1
%
% -------------------------------------------------------------------------
% USGAGE:
% -------------------------------------------------------------------------
% [KH_O2] = gasmolsol(35,20,'O2')
% KH_O2 = 1.1289
%
% Author: David Nicholson dnicholson@whoi.edu
% Also see: gasmolfract.m, gasmoleq.m
% =========================================================================

function [sol] = gasmolsol(SP,pt,gas)

soleq = gasmoleq(SP,pt,gas);
[p_h2o] = vpress(SP,pt);
% water vapour pressure correction
sol = soleq./(gasmolfract(gas).*(1-p_h2o));

