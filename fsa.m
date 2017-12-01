% =========================================================================
% FSA - wrapper function for calculating sea to air gas transfer using a
% specific GE parameterization
%
% [Fd, Fc, Fp, Deq, k] = fas(C_a,C_w,u10,S,T,gas,param,...)
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
%
% [Fd, Fc, Fp, Deq, k] = fsa(C_w,C_a,u10,S,T,gas,param,...)
% T = 10; u10 = 12; S = 35;
% pAr = gasmolfract('Ar')*(1-vpress(35,10));
% USE CONCENTRATION UNITS FOR WATER
% [Fd, Fc, Fp, Deq, k] = fsa(0.014,pAr,u10,S,T,'Ar','L13','w_units','molm3')
% USE DEFAULT UNITS (atm for both w and a)
% [Fd, Fc, Fp, Deq, k] = fsa(380e-6,280e-6,u10,S,T,'CO2','N16')
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% C:        gas concentration (mol m-3)
% u10:      10 m wind speed (m/s)
% S:        Sea surface salinity (PSS)
% T:        Sea surface temperature (deg C)
% gas:      code for gas, formatted as a string, e.g., 'O2' 
% param:    abbreviation for parameterization:
%               Sw07 = Sweeney et al. 2007
%               S09 = Stanley et al. 2009
%               N11 = Nicholson et al. 2011
%               N16 = Nicholson et al. 2016
%               L13 = Liang et al. 2013 
%
%       Code    Gas name        Reference
%       ----   ----------       -----------
%       He      Helium          Weiss 1971
%       Ne      Neon            Hamme and Emerson 2004
%       Ar      Argon           Hamme and Emerson 2004
%       Kr      Krypton         Weiss and Keiser 1978
%       Xe      Xenon           Wood and Caputi 1966
%       N2      Nitrogen        Hamme and Emerson 2004   
%       O2      Oxygen          Garcia and Gordon 1992 
%
% OPTIONAL:
% slp:      sea level pressure (atm) default is 1
% rh:       relative humidity expressed as the fraction of saturation 
%           (0.5 = 50% RH). rh is an optional. default is 1
% -------------------------------------------------------------------------
% OUTPUTS: !!! POSITIVE FLUXES ARE FROM WATER TO AIR !!!
% -------------------------------------------------------------------------
% Fd    Diffusive flux                                (mol m-2 s-1)
% Fc:   Flux from fully collapsing small bubbles      (mol m-2 s-1)
% Fp:   Flux from partially collapsing large bubbles  (mol m-2 s-1)
% Deq:  Equilibrium supersaturation                   (unitless (%sat/100))
% k:    Diffusive gas transfer velocity               (m s-1)
%
% Note: Total air-sea flux is Ft = Fd + Fc + Fp
%
% -------------------------------------------------------------------------
% AUTHOR:
% -------------------------------------------------------------------------
% Author: David Nicholson dnicholson@whoi.edu 
% Version: 29 November 2017
%
% COPYRIGHT:---------------------------------------------------------------
%
% Copyright 2017 Cara Manning and David Nicholson
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
% =========================================================================

function [Fd, Fc, Fp, Deq,k] = fsa(C_w,C_a,u10,S,T,gas,param,varargin)


%% parse inputs

defaultCwUnits = 'atm';
expectedCwUnits = {'atm','molm3'};
defaultCaUnits = 'atm';
expectedCaUnits = {'atm','dryratio'};
defaultSlp = 1;
defaultRh = 1;

%%% parse input parameters
persistent p
if isempty(p)
    p = inputParser;    
    addRequired(p,'C_w',@(x) isnumeric(x));
    addRequired(p,'C_a',@(x) isnumeric(x));
    addRequired(p,'u10',@(x) isnumeric(x));
    addRequired(p,'S',@(x) isnumeric(x));
    addRequired(p,'T',@(x) isnumeric(x));
    addRequired(p,'gas',@isstr);  
    addRequired(p,'param',@isstr);
    addParameter(p,'w_units',defaultCwUnits,@(x) any(validatestring(x,expectedCwUnits)));
    addParameter(p,'a_units',defaultCaUnits,@(x) any(validatestring(x,expectedCaUnits)));
    addParameter(p,'rh',defaultRh,@(x) mean(x) >= 0 & mean(x) <= 1);
    addParameter(p,'slp',defaultSlp,@(x) isnumeric(x));
end
parse(p,C_w,C_a,u10,S,T,gas,param,varargin{:});
inputs = p.Results;


C_w = inputs.C_w;
C_a = inputs.C_a;
u10 = inputs.u10;
S = inputs.S; 
T = inputs.T;
gas = inputs.gas;
param = inputs.param;
w_units = inputs.w_units;
a_units = inputs.a_units;
rh = inputs.rh;
slp = inputs.slp;

if strcmpi(w_units,'atm')
    K0 = gasmolsol(S,T,1,gas);
    C_w_molm3 = K0.*C_w;
elseif strcmpi(w_units,'molm3')
    C_w_molm3 = C_w;
end
if strcmpi(a_units,'atm')
    C_a_atm = C_a;
elseif strcmpi(a_units,'dryratio')
    C_a_atm = C_a.*(slp - rh.*vpress(S,T));
    %C_eq_molm3 = K0.*C_a_atm;
end


switch upper(param)
    case {'W14','SW07','W92A','W92B'}
        [~,Sc] = gasmoldiff(S,T,gas);
        k = kgas(u10,Sc,param);
        Fd = k.*(C_w_molm3 - C_eq_molm3);
        Fc = 0 .* zeros(size(Fd));
        Fp = Fc;
        Deq = Fp;
    case 'N11'
        [Fd, Fc, Fp, Deq, k] = fsa_N11(C_w_molm3,C_a_atm,u10,S,T,gas);
    case 'N16'
        [Fd, Fc, Fp, Deq, k] = fsa_N16(C_w_molm3,C_a_atm,u10,S,T,gas);
    case 'L13'
        [Fd, Fc, Fp, Deq, k] = fsa_L13(C_w_molm3,C_a_atm,u10,S,T,gas);
        
end

% switch upper(param)
%     case 'S09'
%         [Fd, Fc, Fp, Deq, k] = fas_S09(C,u10,S,T,slp,gas,rh);
%     case 'N11'
%         [Fd, Fc, Fp, Deq, k] = fas_N11(C,u10,S,T,slp,gas,rh);
%     case 'N16'
%         [Fd, Fc, Fp, Deq, k] = fas_N16(C,u10,S,T,slp,gas,rh);
%     case 'SW07'
%         [Fd, Fc, Fp, Deq, k] = fas_Sw07(C,u10,S,T,slp,gas,rh);
%     case 'L13'
%         [Fd, Fc, Fp, Deq, k] = fas_L13(C,u10,S,T,slp,gas,rh);
%     case 'I11'
%         [Fd, Fc, Fp, Deq, k] = fas_I11(C,u10,S,T,slp,gas,rh);
%     % TODO: W14
%     otherwise
%         error('Only S09, N11, I11, Sw07 and L13 are supported. See fas_Fd for more parameterizations.');
% end


end