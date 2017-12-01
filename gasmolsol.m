% =========================================================================
% GASMOLSOL.M - calculates Henry's Law solubility for a gas in equilibrium
% with a atmosphere fugacity of 'pgasdry'
% in mol m-3 
%
% This is a wrapper function. See individual solubility functions for more
% details.
%
% [sol] = gasmolsol(SP,pt,pgasdry,gas)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% SP        Practical Salinity
% pt        Potential temperature [degC]
% pgas      gas fugacity [atm] if empty, assumes 1 atm total pressure with
%               saturated water vapor presure
% gas       gas string: He, Ne, Ar, Kr, Xe, O2 or N2
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sol       Henry's Law solubility in mol m-3 atm-1
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% [KH_O2] = gasmolsol(35,20,'O2')
% KH_O2 = 1.1293
%
% Author: David Nicholson dnicholson@whoi.edu
% Also see: gasmolfract.m, gasmoleq.m
% =========================================================================

function [sol] = gasmolsol(SP,pt,pgas,gas)

if isempty(pgas)
    soleq = gasmoleq(SP,pt,gas);
    [p_h2o] = vpress(SP,pt);
    % water vapour pressure correction
    sol = soleq./(gasmolfract(gas).*(1-p_h2o));
else
    sol = 1000.*pgas.*gasBunsen(SP,pt,gas)./gasmolvol(gas);
end

