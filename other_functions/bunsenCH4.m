function [ B, c_umolkg ] = bunsenCH4(S,T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


TK = T + 273.15;
%% For Methane
% REFERENCES:
% Wiesenburg, D. A., and N. L. Guinasso (1979), Equilibrium solubilities of
%   methane, carbon monoxide, and hydrogen in water and sea water, J. Chem. 
%   Eng. Data, 24(4), 356?360, doi:10.1021/je60083a006.
%
% Yamamoto, S., J. B. Alcauskas, and T. E. Crozier (1976), Solubility of 
%   methane in distilled water and seawater, J. Chem. Eng. Data, 21(1), 
%   78?80, doi:10.1021/je60068a029.


switch gas
% These constants are refit to Yamamoto et al. 1976 found in Wiesenburg and
% Guinasso 1979
    case {'CH4','methane'}
        A1 = -68.8862;
        A2 = 101.4956;
        A3 = 28.7314;
        B1 = -0.076146;
        B2 = 0.043970;
        B3 = -0.0068672;       
    case {'CO','carbon monoxide'}
        A1 = -47.6148;
        A2 = 69.5068;
        A3 = 18.7397;
        B1 = 0.045657;
        B2 = -0.040721;
        B3 = 0.0079700;
    case {'H2','hydrogen'}
        A1 = -47.8948;
        A2 = 65.0368;
        A3 = 20.1709;
        B1 = -0.082225;
        B2 = 0.049564;
        B3 = -0.0078689;
    case {'CH4yamamato'} % original constants from Yamamoto et al.
        A1 = -67.1952;
        A2 = 99.1624;
        A3 = 27.9015;
        B1 = -0.072909;
        B2 = 0.041674;
        B3 = -0.0064603;
end
        
logB = A1 + A2.*(100./TK) + A3.*log(TK./100) + S.*( B1 + B2.*(T./100) ...
    + B3.*((T./100).^2));

B = exp(logB);


% Eqn (7) of Wiesenburg 1979

% % A1 = -417.5053;
% % A2 =  599.8626;
% % A3 =  380.3636;
% % A4 = -62.0764;  %This term corrects for the use of dry air
% % B1 = -0.064236;
% % B2 = 0.034980;
% % B3 = -0.0052732;
% %lnC = log(pGdry) + A1 + A2.*(100./TK) + A3.*log(TK./100) + A4.*(TK./100)...
% %    + S.*( B1 + B2.*(T./100) + B3.*((T./100).^2));
% % c_umolkg = exp(lnC)./1000;


end

