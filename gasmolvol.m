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
%
% AUTHOR:------------------------------------------------------------------
% David Nicholson dnicholson@whoi.edu
% Woods Hole Oceanographic Institution
%==========================================================================

function vmol = gasmolvol(gas)
gasn = lower(gas);
switch gasn
    case 'he'
        vmol = 22.426;
    case 'ne'
        vmol = 22.425;
    case 'ar'
        vmol = 22.393;
    case 'kr'
        vmol = 22.352;
    case 'xe'
        vmol = 22.258;
    case 'o2'
        vmol = 22.392;
    case 'n2'
        vmol = 22.404;
    case 'co2'
        % from http://cdiac.ornl.gov/ftp/cdiac74/sop24.pdf
        vmol = 22.414.*0.99498;
    case 'ch4'
        vmol = 22.360;
    case 'h2'
        vmol = 22.428;
    case 'n2o'
        vmol = 22.243;
    otherwise
        error('Gas name must be He, Ne, Ar, Kr, Xe, N2, O2, CO2, CH4, H2 or N2O');
end
end               