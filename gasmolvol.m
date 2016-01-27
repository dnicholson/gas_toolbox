% GASMOLVOL.M - Function giving molar volumes at standard temperature and
% pressure (STP), which is defined as 273.15 K, 1 atm of the pure gas
%
% USAGE:-------------------------------------------------------------------
% vmol = gasmolvol(gas)
%
% INPUT:-------------------------------------------------------------------
% gas       gas string: He, Ne, Ar, Kr, Xe, O2, or N2
%
% REFERENCE:---------------------------------------------------------------
% Dymond, J. H., Marsh, K. N., Wilhoit, R. C., & Wong, K. C. (2002). Virial
% Coefficients of Pure Gases and Mixtures. Landolt-Börnstein Numerical Data
% and Functional Relationships in Science and Technology New Series Group
% IV. Physical Chemistry Vol, 21.
%
% AUTHOR:------------------------------------------------------------------
% David Nicholson dnicholson@whoi.edu
% Woods Hole Oceanographic Institution
%==========================================================================

function vmol = gasmolvol(gas)
if strcmpi(gas,'He')
    vmol = 22.426;
elseif strcmpi(gas,'Ne')
    vmol = 22.425;
elseif strcmpi(gas,'Ar')
    vmol = 22.393;
elseif strcmpi(gas,'Kr')
    vmol = 22.352;
elseif strcmpi(gas,'Xe')
    vmol = 22.258;
elseif strcmpi(gas,'O2')
    vmol = 22.392;
elseif strcmpi(gas,'N2')
    vmol = 22.404;
elseif strcmpi(gas,'CO2')
    % from http://cdiac.ornl.gov/ftp/cdiac74/sop24.pdf
    vmol = 22.414.*0.99498;
else
    error('Gas name must be Ne, Ar, Kr, Xe, N2, or O2');
end               