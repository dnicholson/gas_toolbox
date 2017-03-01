% =========================================================================
% [sol] = gasmoleq(SP,pt,gas)
%
% GASMOLEQ.M - calculates equilibrium solubility of a dissolved gas
% in mol/m^3 at an absolute pressure of 101325 Pa (sea pressure of 0 
% dbar) including saturated water vapor. 
%
% This is a wrapper function. See individual solubility functions for more
% details.
%
% This function uses the GSW Toolbox solubility functions when available,
% except for Ne. There was a bug in gsw_Nesol_SP_pt and gsw_Nesol, in
% versions downloaded prior to Sept 23, 2015. Therefore, we provide a
% correct solubility function to avoid errors for people with older
% versions of the GSW Toolbox.
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% [O2eq]    = gasmoleq(35,20,'O2')
% O2eq      = 0.2311
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% SP        Practical salinity [PSS-78]
% pt        Potential temperature [degC]
% gas       gas string: He, Ne, Ar, Kr, Xe, O2 or N2
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sol       gas equilibrium solubility in mol/m^3
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  See also the references within each solubility function. 
%
% -------------------------------------------------------------------------
% AUTHORS:
% -------------------------------------------------------------------------
% Cara Manning, cmanning@whoi.edu
% David Nicholson, dnicholson@whoi.edu
%
% -------------------------------------------------------------------------
% LICENSE:
% -------------------------------------------------------------------------
% Copyright 2015 Cara Manning and David Nicholson 
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
% =========================================================================

function [sol] = gasmoleq(SP,pt,gas)

% Check for GSW solubility functions
% Commented out to speed up performance
% if exist('gsw_Arsol_SP_pt','file') == 0
%     error('gsw_Arsol_SP_pt not in MATLAB path. Download GSW Toolbox at http://www.teos-10.org/');
% end;

% Calculate potential density at surface
SA = SP.*35.16504./35;
CT = gsw_CT_from_pt(SA,pt);
rho = gsw_sigma0(SA,CT)+1000;

gasn = lower(gas);
% calculate equilibrium solubility gas concentration in micro-mol/kg
if strcmpi(gasn, 'he')
    sol_umolkg = gsw_Hesol_SP_pt(SP,pt);
elseif strcmpi(gasn, 'ne')
    % bug in gsw_Nesol... v. 3.05 prior Sept 23, 2015 
    % the bug returned solubility in nmol kg-1 instead of umol kg-1
    sol_umolkg = Nesol(SP,pt);
elseif strcmpi(gasn, 'ar')
    sol_umolkg = gsw_Arsol_SP_pt(SP,pt);
elseif strcmpi(gasn, 'kr')
    sol_umolkg = gsw_Krsol_SP_pt(SP,pt);
elseif strcmpi(gasn, 'xe')
    sol_umolkg = Xesol(SP,pt);
elseif strcmpi(gasn, 'n2')
    sol_umolkg = gsw_N2sol_SP_pt(SP,pt);
elseif strcmpi(gasn, 'o2')
    sol_umolkg = gsw_O2sol_SP_pt(SP,pt);
elseif strcmpi(gasn, 'n2o')
    sol_umolkg = gsw_N2Osol_SP_pt(SP,pt);
else
    error('Gas name must be He, Ne, Ar, Kr, Xe, O2, N2O or N2');
end

% convert from micro-mol/kg to mol/m3
sol = rho.*sol_umolkg./1e6;

