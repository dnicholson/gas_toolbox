function [kw_wt] = kw_weighting_ice(kw,dt,ndays,zml,ice)

%-------------------------------------------------------------------------
% ABOUT:
% Calculate the weighted piston velocity for historic wind data; applies
% ice coverage correction
% Time weighting following Teeter et al. (2018) equation 5 (modification to
% Reuer et al. (2007)
% 
% INPUT:
% kw = vector of piston velocities [m/day], from a single location. NOTE:
% time is forward -- kw(1) is the oldest piston velocity; kw(end) is "today"
% dt = time interval [day]
% ndays = total number of days overwhich to weight the piston velocity
% zml = mixed layer depth (vector, same size as kw), or scalar
% ice = ice coverage % (vector, same size as kw), or scalar
% 
% OUTPUT:
% kw_wt = weighted kw [m/day]
% 
% Script created by:
% R. Izett 
% Last updated: Apr 2020 by Cara Manning
% UBC Oceanography
%
% REFERENCES

%-------------------------------------------------------------------------

%--- Flip kw, zml and ice vectors 
     % such that X(1) is the most recent observations
    kw = flip(kw);
    zml = flip(zml);
    ice = flip(ice);

%--- Apply ice correction to each k value
    kw = kw .* (1-ice); %as per Butterworth & Miller, 2016

%--- Create weighting and fraction ML ventilated vectors
    wt = nan(size(kw)); %weighting
    fr = nan(size(kw)); %fraction ventilated
    
    wt(1) = 1; %assign weighting of 1 to most recent observation
    fr(1) = kw(1).*dt./zml(1); %calculate fraction ventilated on day of sampling

%--- go backwards and fill in vectors
    for ff = 2:ndays/dt+1
        fr(ff) = kw(ff).*dt./zml(ff);
        if fr(ff)>=1; fr(ff)=.9999; end
        wt(ff) = wt(ff-1).*(1-fr(ff-1));
    end
  
%--- Calculate weighted kw
    prod = kw.*wt; %the product of kw * weight for each observation
    
    kw_wt = (sum(prod)) / (sum(wt));
    
end
    
    
    


