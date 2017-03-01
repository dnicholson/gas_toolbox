% =========================================================================
% KGAS - gas transfer coefficient for a range of windspeed-based
% parameterizations
%
% [kv] = kgas(u10,Sc,param)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% u10       10-m wind speed (m/s)
% Sc        Schmidt number 
% param     abbreviation for parameterization:
%           W14  = Wanninkhof 2014
%           W92a = Wanninkhof 1992 - averaged winds
%           W92b = Wanninkhof 1992 - instantaneous or steady winds
%           Sw07 = Sweeney et al. 2007
%           Ho06 = Ho et al. 2006
%           Ng00 = Nightingale et al. 2000
%           LM86 = Liss and Merlivat 1986
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% kv       Gas transfer velocity in m s-1
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% k = kgas(10,1000,'W14')
% k = 5.6643e-05
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
% Wanninkhof, R. (2014). Relationship between wind speed and gas exchange
% over the ocean revisited. Limnol. Oceanogr. Methods, 12(6), 351-362.
%
% Wanninkhof, R. (1992). Relationship between wind speed and gas exchange
% over the ocean. J. Geophys. Res, 97(25), 7373-7382.
%
% Sweeney, C., Gloor, E., Jacobson, A. R., Key, R. M., McKinley, G.,
% Sarmiento, J. L., & Wanninkhof, R. (2007). Constraining global air-sea
% gas exchange for CO2 with recent bomb 14C measurements. Global
% Biogeochem. Cy.,21(2).
%
% Ho, D. T., Law, C. S., Smith, M. J., Schlosser, P., Harvey, M., & Hill,
% P. (2006). Measurements of air?sea gas exchange at high wind speeds in
% the Southern Ocean: Implications for global parameterizations. Geophys.
% Res. Lett., 33(16).
%
% Nightingale, P. D., Malin, G., Law, C. S., Watson, A. J., Liss, P. S.,
% Liddicoat, M. I., et al. (2000). In situ evaluation of air-sea gas
% exchange parameterizations using novel conservative and volatile tracers.
% Global Biogeochem. Cy., 14(1), 373-387.
%
% Liss, P. S., & Merlivat, L. (1986). Air-sea gas exchange rates:
% Introduction and synthesis. In The role of air-sea exchange in
% geochemical cycling (pp. 113-127). Springer Netherlands.
%
% -------------------------------------------------------------------------
% AUTHORS
% -------------------------------------------------------------------------
% David Nicholson dnicholson@whoi.edu, Woods Hole Oceanographic Institution
% Modified by Cara Manning cmanning@whoi.edu
% Version 2.0
%
% -------------------------------------------------------------------------
% COPYRIGHT
% -------------------------------------------------------------------------
% Copyright 2015 David Nicholson and Cara Manning 
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

function [kv] = kgas(u10,Sc,param)

quadratics = {'W14','W92a','W92b','Sw07','Ho06'};
% should be case insensitive
parup = upper(param);
if ismember(parup,upper(quadratics))
    if strcmpi(parup,'W14')
        A = 0.251;
    elseif strcmpi(parup,'W92A')
        A = 0.39;
    elseif strcmpi(parup,'W92B')
        A = 0.31;
    elseif strcmpi(parup,'SW07')
        A = 0.27;
    elseif strcmpi(parup,'HO06')
        A = 0.254; %k_600 = 0.266
    end
    k_cm = A*u10.^2.*(Sc./660).^-0.5; %cm/h
    kv = k_cm./(100*60*60); %m/s
elseif strcmpi(parup,'NG00')
    k600 = 0.222.*u10.^2 + 0.333.*u10;
    k_cm = k600.*(Sc./600).^-0.5; % cm/h
    kv = k_cm./(100*60*60); % m/s
elseif strcmpi(parup,'LM86')
    k600 = zeros(1,length(u10));
    
    l = find(u10 <= 3.6);
    k600(l) = 0.17.*u10(l);
    
    m = find(u10 > 3.6 & u10 <= 13);
    k600(m) = 2.85.*u10(m)-9.65;
    
    h = find(u10 > 13);
    k600(h) = 5.9.*u10(h)-49.3;
    
    k_cm = k600.*(Sc./600).^-0.5;
    k_cm(l) = k600(l).*(Sc./600).^(-2/3);
    kv = k_cm./(100*60*60); % m/s
else
    error('parameterization not found');
end
end