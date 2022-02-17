function [kw_wt,wt] = kw_weighting(kw,dt,ndays,zml)

%-------------------------------------------------------------------------
% ABOUT:
% Calculate the weighted piston velocity for historic wind data.
% 
% INPUT:
% kw = vector of piston velocities [m/day], from a single location. NOTE:
% time is forward -- kw(1) is the oldest piston velocity; kw(end) is "today"
% dt = time interval [day]
% ndays = total number of days overwhich to weight the piston velocity
% zml = mixed layer depth (vector, same size as kw), or scalar
% 
% OUTPUT:
% kw_wt = weighted kw [m/day]
% wt = weighting 
% 
% NOTE: Kw weighting based on Reuer et al., 2007, with corrections from
% Teeter et al. (2018).
% 
% Script created by:
% R. Izett 
% Last updated: Apr 2020 by Cara Manning
% Adapted from time_weigted_kw_kaiser by M. Long (2007).
%-------------------------------------------------------------------------

%--- Flip kw and zml vectors 
     % such that kw(1) and zml(1) are the most recent observations
    kw = flip(kw);
    zml = flip(zml);

%--- Create weighting and fraction ML ventilated vectors
    wt = nan(size(kw)); %weighting
    fr = nan(size(kw)); %fraction ventilated
    
    wt(1) = 1; %assign weighting of 1 to "today's"/most recent observation
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
