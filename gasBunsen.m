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
% See references for individual gas solubility functions.
%
% AUTHORS:-----------------------------------------------------------------
% Cara Manning (cmanning@whoi.edu) Woods Hole Oceanographic Institution
% David Nicholson (dnicholson@whoi.edu)
% Version: 1.0 // September 2015
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2015 Cara Manning 
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beta=gasBunsen(SP,pt,gas)

pdry = 1 - vpress(SP,pt); % pressure of dry air for 1 atm total pressure

% equilib solubility in mol/m3
Geq = gasmoleq(SP,pt,gas);

% calc beta 
beta = Geq.*(gasmolvol(gas)./1000)./(pdry.*gasmolfract(gas));
end
