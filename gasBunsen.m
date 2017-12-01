% beta = gasBunsen(SP,pt,gas)
% Function to calculate Bunsen coefficient
%
% USAGE:-------------------------------------------------------------------
% beta=gasBunsen(SP,pt,gas)
%
% DESCRIPTION:-------------------------------------------------------------
% Calculate the Bunsen coefficient, which is defined as the volume of pure
% gas at standard temperature and pressure (0 degC, 1 atm) that will
% dissolve into a volume of water at equilibrium when exposed to a
% pressure of 1 atm of the gas. 
%
% Side note:  IUPAC STP is defined as 0 deg C, 1 bar rather than 1 atm. 
% (1atm = 1.01325 bar).  However, prior to ~1982 1 atm was commonly used to
% define STP.  Care must be taken to identify which pressure unit is used
% in "STP" for Bunsen coefficient.  Traditionally (e.g., Weiss) 1 atm is  
% used for Bunsen, not 1 bar. For consistency with literature, we use 1 atm 
% for "STP"
%
%
% INPUTS:------------------------------------------------------------------
% SP:    Practical salinity (PSS)
% pt:    Potential temperature (deg C)
% gas:  code for gas (He, Ne, Ar, Kr, Xe, N2, or O2)
%
% OUTPUTS:-----------------------------------------------------------------
% beta: Bunsen coefficient                  (L gas)/(L soln * atm gas)
%
% REFERENCE:---------------------------------------------------------------
%
% % REFERENCES:
% CH4,CO and H2:
% Wiesenburg, D. A., and N. L. Guinasso (1979), Equilibrium solubilities of
%   methane, carbon monoxide, and hydrogen in water and sea water, J. Chem. 
%   Eng. Data, 24(4), 356?360, doi:10.1021/je60083a006.
% CH4:
% Yamamoto, S., J. B. Alcauskas, and T. E. Crozier (1976), Solubility of 
%   methane in distilled water and seawater, J. Chem. Eng. Data, 21(1), 
%   78?80, doi:10.1021/je60068a029.
%
% N2O:
% Weiss, R. F., and B. A. Price (1980), Nitrous oxide solubility in water 
%   and seawater, Mar. Chem., 8(4), 347?359, doi:10.1016/0304-4203(80)90024-9.
%
% CO2:
% Weiss, R. F. 1974. Carbon dioxide in water and seawater: the solu-
% bility of a non-ideal gas. Mar. Chem. 2:203-215 
% doi:10.1016/0304-4203(74)90015-2. (in Wanninkhof 2014)
%
% He,Ne,Ar,Kr,Xe,N2,O2:
% See references for individual gas solubility functions.
%
%
% AUTHORS:-----------------------------------------------------------------
% Cara Manning (cmanning@whoi.edu) Woods Hole Oceanographic Institution
% David Nicholson (dnicholson@whoi.edu)
% Version: 1.0 // September 2015
% Version: 2.0 // December 2016
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2015-2016 David Nicholson, Cara Manning 
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beta=gasBunsen(SP,pt,gas)

gasn = lower(gas);
if ismember(gasn,{'he','ne','ar','kr','xe','n2','o2'})
pdry = 1 - vpress(SP,pt); % pressure of dry air for 1 atm total pressure

% equilib solubility in mol/m3
Geq = gasmoleq(SP,pt,gas);

% calc beta 
beta = Geq.*(gasmolvol(gas)./1000)./(pdry.*gasmolfract(gas));

% The following gases use fit coefs for beta directly
elseif ismember(gasn,{'ch4','co','h2','ch4yamamato','co2','n2o'})
    ufac = 1;  % unit conversion if needed
    switch gasn
        % These constants are refit to Yamamoto et al. 1976 found in Wiesenburg and
        % Guinasso 1979
        case {'ch4','methane'}
            A1 = -68.8862;
            A2 = 101.4956;
            A3 = 28.7314;
            B1 = -0.076146;
            B2 = 0.043970;
            B3 = -0.0068672;
        case {'co','carbon monoxide'}
            A1 = -47.6148;
            A2 = 69.5068;
            A3 = 18.7397;
            B1 = 0.045657;
            B2 = -0.040721;
            B3 = 0.0079700;
        case {'h2','hydrogen'}
            A1 = -47.8948;
            A2 = 65.0368;
            A3 = 20.1709;
            B1 = -0.082225;
            B2 = 0.049564;
            B3 = -0.0078689;
        case {'n2o','nitrous oxide'}
            A1 = -62.7062;
            A2 = 97.3066;
            A3 = 24.1406;
            B1 = -0.058420;
            B2 = 0.033193;
            B3 = -0.0051313;
            pt = pt.*1.00024; % potential temperature in degress C on 
              % the 1968 International Practical Temperature Scale IPTS-68.
            ufac = gasmolvol('N2O');  % convert mol/L --> L(stp)/L(sol)
        case {'ch4yamamato'}
            A1 = -67.1952;
            A2 = 99.1624;
            A3 = 27.9015;
            B1 = -0.072909;
            B2 = 0.041674;
            B3 = -0.0064603;
        case {'co2'}
            A1 = -58.08931;
            A2 = 90.5069;
            A3 = 22.2940;
            B1 = 0.027766;
            B2 = -0.025888;
            B3 = 0.0050578;
            ufac = gasmolvol('CO2'); % convert from mol L-1 atm-1
    end
    TK = pt + gsw_T0;
    logB = A1 + A2.*(100./TK) + A3.*log(TK./100) + SP.*( B1 + B2.*(TK./100) ...
        + B3.*((TK./100).^2));
    
    beta = ufac.*exp(logB);
end
end
