% u10 = calc_u10(umeas,hmeas)
%
% USAGE:-------------------------------------------------------------------
% 
% [u10] = calc_u10(5,4)
%
% >u10 = 5.5302
%
% DESCRIPTION:-------------------------------------------------------------
% Scale wind speed from measurement height to 10 m height
%
% INPUTS:------------------------------------------------------------------
% umeas:  measured wind speed (m/s)
% hmeas:  height of measurement (m), can have dimension 1x1 or same as umeas
%
% OUTPUTS:-----------------------------------------------------------------
%
% Fs:   Surface gas flux                              (mol m-2 s-1)
% Fp:   Flux from partially collapsing large bubbles  (mol m-2 s-1)
% Fc:   Flux from fully collapsing small bubbles      (mol m-2 s-1)
% Deq:  Equilibrium supersaturation                   (unitless (%sat/100))
%
% REFERENCE:---------------------------------------------------------------
%
% Hsu S, Meindl E A and Gilhousen D B (1994) Determining the Power-Law
%    Wind-Profile Exponent under Near-Neutral Stability Conditions at Sea 
%    J. Appl. Meteor. 33 757-765
%
% AUTHOR:---------------------------------------------------------------
% Cara Manning cmanning@whoi.edu
% MIT/WHOI Joint Program in Oceanography
% Version: 1.0 // September 2015
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2015 Cara Manning 
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u10 = calc_u10(umeas,hmeas)

u10 = umeas.*(10./hmeas).^0.11;
end