% -------------------------------------------------------------------------
% DEMO_WEIGHTED_GAS_FLUXES.M
% Calculate sea-air flux of CH4 and N2O by three different methods:
% 1) instantaneous flux
% 2) 30-day weighted flux
% 3) 60-day weighted flux
% All methods use the Ho et al. (2006) parameterization for gas exchange
% and assume that the gas transfer velocity scales linearly as a function
% of the fraction of open water (Butterworth and Miller, 2016)
% The weighting is performed using the Teeter et al. (2018) method which
% modifies equations first derived in Reuer et al. (2007)
% Positive flux is from the sea to the air.
%
% Author: Cara Manning, adapted from code by Robert Izett
%
% -------------------------------------------------------------------------
% REFERENCES
% Butterworth, B. J., and Miller, S. D. (2016). Air‐sea exchange of carbon
% dioxide in the Southern Ocean and Antarctic marginal ice zone.
% Geophysical Research Letters, 43(13), 7223-7230.
% https://doi.org/10.1002/2016GL069581
% 
% Ho, D. T., Law, C. S., Smith, M. J., Schlosser, P., Harvey, M., and Hill,
% P. (2006), Measurements of air-sea gas exchange at high wind speeds in
% the Southern Ocean: Implications for global parameterizations, Geophys.
% Res. Lett., 33, L16611, doi:10.1029/2006GL026817.
%
% Reuer, M. K., Barnett, B. A., Bender, M. L., Falkowski, P. G., and
% Hendricks, M. B. (2007). New estimates of Southern Ocean biological
% production rates from O2/Ar ratios and the triple isotope composition of
% O2. Deep Sea Research Part I: Oceanographic Research Papers, 54(6),
% 951-974. https://doi.org/10.1016/j.dsr.2007.02.007
%
% Teeter, L., Hamme, R. C., Ianson, D., & Bianucci, L. (2018). Accurate
% estimation of net community production from O2/Ar measurements. Global
% Biogeochemical Cycles, 32, 1163– 1181.
% https://doi.org/10.1029/2017GB005874
%
% -------------------------------------------------------------------------
% STEPS TO CALCULATE WEIGHTED GAS FLUXES
% 1) Following DEMO_DOWNLOAD.M, download gridded wind, sea level pressure,
% and sea ice data covering the lat/lon range of all samples and going back
% in time at least 30 or 60 days prior to each observation. The data for
% each variable should be saved in an array structure (struct).
% The struct field names will be mdate, lat, lon and the variable name
% (slp, wind, or ice_frac)
% The field with the variable data will have units time x lat x lon
%
% If you are working with an area with no ice cover, you can set ice_frac
% to 0 or modify the script to use kw_weighting instead of
% kw_weighting_ice.
%
% 2) Prepare an array structure containing the gas data for each sample.
% A demo gas data file DEMO_GAS_DATA.CSV is included in the DATA folder for
% testing purposes.
% For CH4 and N2O fluxes as calculated here, the gas data file must
% include:
% gd.lat = latitude
% gd.lon = longitude
% gd.temp = temperature (deg C)
% gd.sal = practical salinity (PSS-78)
% gd.press = pressure (db)
% gd.depth = depth (m)
% gd.ch4_nmolkg = CH4 concentration in nmol/kg
% gd.ch4_eq_nmolkg = CH4 equilibrium concentration in nmol/kg, calculated
% using CH4sol.m
% gd.n2o_nmolkg = N2O concentration in nmol/kg
% gd.n2o_eq_nmolkg = N2O equilibrium concentration in nmol/kg, calculated
% using N2Osol.m
% gd.mld = mixed layer depth (m), calculated using calcmld.m or similar

% 3) Run the following code to calculate the gas fluxes for each sample and
% save to an array structure. The values are also plotted.
%
% The function distance.m is part of the Matlab mapping toolbox. If you do
% not have access to this toolbox, you can modify the code to use
% gsw_distance (from the GSW Toolbox) or m_lldist (from M_Map, available at
% https://www.eoas.ubc.ca/~rich/map.html).
% -------------------------------------------------------------------------
clc; clear; close all;

% load in the data downloaded through demo_download.m
% change ndir to match your directory structure
ndir = 'C:\Users\cmanning\Dropbox\gas_toolbox\weighted_flux\data\';

load([ndir,'demo_SLP.mat']);
load([ndir,'demo_CCMP.mat']);
load([ndir,'demo_AMSR2.mat']);

demo_gas_data = readtable('demo_gas_data.csv');

%%
wt_t_30 = 30; % weighting time of 30 days
wt_t_60 = 60; % weighting time of 60 days


mbar_to_atm = 1./1013.25; % conversion factor for pressure data
                          % NCEP provides pressure in mbar but solubility
                          % functions use atm


% rename the imported files
satice = demo_AMSR2; 
satwind = demo_CCMP;
satslp = demo_SLP;
gd = demo_gas_data;

% convert slp from mbar to atm
satslp.slp = satslp.slp .* mbar_to_atm; % slp in atm

% save the gas data to new variable names
s_time = datenum(gd.yyyy,gd.mm,gd.dd,gd.HH,gd.MM,gd.SS);
s_lat = gd.lat;
s_lon = gd.lon;
s_mld = gd.mld;
s_T = gd.temp;
s_S = gd.sal;
s_P = gd.press;
s_ch4 = gd.ch4_nmolkg;
s_n2o = gd.n2o_nmolkg;
s_ch4_eq = gd.ch4_eq_nmolkg;
s_n2o_eq = gd.n2o_eq_nmolkg;
s_station = gd.station;
s_depth = gd.depth;

% create structure to save the interpolated data
und.time=s_time; und.lat = s_lat; und.lon= s_lon;

%--- Wind matrix
% Interpolate 3D wind matrix to the observations
[X_wind,Y_wind,Z_wind] = ndgrid(satwind.mdate,satwind.lat,satwind.lon);
V_wind = satwind.wind;
% calculate instantaneous wind speeds for time of each sample
und.wind = interpn(X_wind,Y_wind,Z_wind,V_wind,und.time,und.lat,und.lon);

%--- SLP matrix
% Interpolate 3D SLP matrix to the observations
[X_slp,Y_slp,Z_slp] = ndgrid(satslp.mdate,satslp.lat,satslp.lon);
V_slp = satslp.slp;
und.slp = interpn(X_slp,Y_slp,Z_slp,V_slp,und.time,und.lat,und.lon); %instantaneous slp for time of cruise

%-------------- instantaneous and 30-day weighted fluxes -----------------
%--- Create historical wind and slp matrices for the samples using the
% weighting period wt_t_30 (30 days) 
wt_t = wt_t_30;
int = 6; % interval = 6 hr for wind and slp
lag = 1:(wt_t *(24/int)); % # observations back in time (1 : # days * # observations per day)
windmat = nan(length(und.time),length(lag)+1); % empty matrix; note: the +1 is to hold the current (i.e. 0 hr lag) observation
t_back = nan(length(und.time),length(lag)+1);

slp_t = wt_t_30;
lag_slp = 1:(slp_t *(24/int));
t_back_slp = nan(length(und.time),length(lag)+1);
slpmat = nan(length(und.time),length(lag_slp)+1); %empty matrix; note: the +1 is to hold the current (i.e. 0 hr lag) observation


for kk=fliplr(lag)
    t_back(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
    wind_back = interpn(X_wind,Y_wind,Z_wind,V_wind,t_back(:,kk),und.lat,und.lon); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid
    windmat(:,lag(end)-lag(kk)+1) = wind_back; %fill matrix
end
t_back(:,end) = und.time;
und.windmat = windmat;
windmat(:,end) = und.wind; %add instantaneous wind

for kk=fliplr(lag_slp)
    t_back_slp(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
    slp_back = interpn(X_slp,Y_slp,Z_slp,V_slp,t_back_slp(:,kk),und.lat,und.lon); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid
    slpmat(:,lag(end)-lag(kk)+1) = slp_back; %fill matrix
end

slpmat(:,end) = und.slp; %add instantaneous slp
und.slpmat = slpmat;

% clear variables no longer needed
clear  X_wind Y_wind Z_wind V_wind ...
    X_slp Y_slp Z_slp V_slp ...
    kk t_back wind_back uw_wind

% --- ice_frac matrix
% for every sample, find the index of the nearest pixel in the ice data
ice_row = nan.*s_lat;
ice_col = nan.*s_lat;
for m = 1:numel(s_lat)
    [arclen,~] = distance(s_lat(m),s_lon(m),satice.lat,satice.lon);
    [~,I] = min(arclen(:));
    [ice_row(m), ice_col(m)] = ind2sub(size(arclen),I);
end

% for every sample, calculate the instantaneous ice conc by
% interpolating to the observation time
und.ice_frac = nan.*und.lat;
for m = 1:numel(ice_row)
    und.ice_frac(m) = interpn(satice.mdate,satice.ice_frac(:,ice_row(m),ice_col(m)),und.time(m));
end

% for any nan values, set ice_frac to 0
und.ice_frac(isnan(und.ice_frac)) = 0;

%--- Create historical ice_frac matrix for the samples using the
% weighting period wt_t_30 (30 days) 
% int = 6; %6 hr interval to match CCMP winds
% lag = 1:(wt_t *(24/int)); % # observations back in time (1 : # days * # observations per day)
icemat = nan(length(und.time),length(lag)+1); % empty matrix; note: the +1 is to hold the current (i.e., 0 hr lag) observation
t_back = nan(length(und.time),length(lag)+1);

ice_back = nan(numel(ice_row),numel(lag)+1);
for kk=fliplr(lag)
    t_back(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g., if kk = 120 and int = 6 hr ==> lag backwards 720 hr (30 d)
    for m = 1:numel(ice_row)
        ice_back(m,kk) = interpn(satice.mdate,satice.ice_frac(:,ice_row(m),ice_col(m)),t_back(m,kk));
    end
    icemat(:,lag(end)-lag(kk)+1) = ice_back(:,kk); %fill matrix
end

icemat(:,end) = und.ice_frac;
und.icemat = icemat;
und.icemat(isnan(und.icemat)) = 0; % for any nan values within the matrix, set to 0 ice cover;

clear icemat ice_back t_back; % clear variables no longer needed

%--- make matrices of T, S, MLD, and Schmidt number
% Assume constant backwards in time
Tmat = repmat(s_T,1,length(lag)+1);
Smat = repmat(s_S,1,length(lag)+1);
zMLmat = repmat(s_mld,1,length(lag)+1);

[~,Scmat_CH4] = gasmoldiff(Smat,Tmat,'CH4');
[Scmat_N2O] = calc_schmidt(Smat,Tmat,'n2o');

% calculate open water gas transfer velocity following Ho et al. (2006)
spd = 60*60*24; % seconds per day
k.CH4 = kgas(windmat,Scmat_CH4,'Ho06').* spd; % m/d
k.N2O = kgas(windmat,Scmat_N2O,'Ho06') .* spd; % m/d

ice = und.icemat;

slp_inst = slpmat(:,end); % instantaneous slp

% calculate instantaneous gas transfer velocity corrected for ice cover
gd.k_inst_CH4 = k.CH4(:,1).*(1-ice(:,1)); % instantaneous CH4 gas tranfer velocity in m/d
gd.k_inst_N2O = k.N2O(:,1).*(1-ice(:,1)); % instantaneous N2O gas tranfer velocity in m/d

% instantenous sea-air gas flux corrected for ice cover and slp
gd.F_CH4_inst = gd.k_inst_CH4.*(s_ch4-s_ch4_eq.*slp_inst)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_inst = gd.k_inst_N2O.*(s_n2o-s_n2o_eq.*slp_inst)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d

% calculate 30-day weighted gas transfer velocity corrected for ice cover
gd.k_wt_30_CH4 = nan.*s_time; % initialize variable
gd.k_wt_30_N2O = nan.*s_time; % initialize variable

for kk = 1:length(Scmat_CH4(:,1))
    gd.k_wt_30_CH4(kk) = kw_weighting_ice(k.CH4(kk,:), int/24, wt_t, zMLmat(kk,:),ice(kk,:)./100); % m/d
    gd.k_wt_30_N2O(kk) = kw_weighting_ice(k.N2O(kk,:), int/24, wt_t, zMLmat(kk,:),ice(kk,:)./100); % m/d
end

clear kk

slp_30 = mean(slpmat,2); % 30-day average slp

ice_inst = ice(:,end); % instantaneous ice fraction
ice_30 = mean(ice,2); % 30-day average ice fraction

% calculate 30-day weighted sea-air flux of gases
% Units are umol/m2/d = m/d * umol/kg * kg/m3
gd.F_CH4_30 = gd.k_wt_30_CH4.*(s_ch4-s_ch4_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_30 = gd.k_wt_30_N2O.*(s_n2o-s_n2o_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
disp('processed 30-day weighted gas fluxes')
   

%-------------- 60-day weighted fluxes -----------------
%--- calculate fluxes using weighting period wt_t_60
%--- Wind matrix
% Fit 3D wind matrix to a 3D cruise track line
[X,Y,Z]=ndgrid(satwind.mdate,satwind.lat,satwind.lon); %make grid from wind data
V = satwind.wind;
und.time=s_time; und.lat = s_lat; und.lon= s_lon;
und.wind = interpn(X,Y,Z,V,und.time,und.lat,und.lon); %instantaneous wind speeds for time of cruise

%Now make historic wind matrix along cruise track (for estimating weighted piston velocity)
wt_t = wt_t_60;
int = 6; %6 hr interval for CCMP and NCEP winds -- if you have daily measurements of wind speed, this would be 24
lag = 1:(wt_t *(24/int)); % # observations back in time (1: # days * # observations per day)
windmat = nan(length(und.time),length(lag)+1); %empty matrix; note: the +1 is to hold the current (i.e. 0 hr lag) observation
t_back = nan(length(und.time),length(lag)+1);

for kk=fliplr(lag)
    t_back(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
    wind_back = interpn(X,Y,Z,V,t_back(:,kk),und.lat,und.lon); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid
    windmat(:,lag(end)-lag(kk)+1) = wind_back; %fill matrix
end
windmat(:,end) = und.wind; %add instantaneous wind
t_back(:,end) = und.time;
und.windmat = windmat;

clear  X Y Z kk t_back wind_back uw_wind

% skip the step to find the nearest pixel as this was already done for the
% 30-day weighting
% % for every sample, find the index of the nearest pixel in the ice data
% ice_row = nan.*s_lat;
% ice_col = nan.*s_lat;
% for m = 1:numel(s_lat)
%     [arclen,~] = distance(s_lat(m),s_lon(m),satice.lat,satice.lon);
%     [~,I] = min(arclen(:));
%     [ice_row(m), ice_col(m)] = ind2sub(size(arclen),I);
% end
% 
% % for every sample, calculate the instantaneous ice conc
% und.ice_frac = nan.*und.lat;
% for m = 1:numel(ice_row)
%     und.ice_frac(m) = interpn(satice.mdate,satice.ice_frac(:,ice_row(m),ice_col(m)),und.time(m));
% end
% 
% % for any nan values, set ice cover to 0%
% und.ice_frac(isnan(und.ice_frac)) = 0;

%--- Create historical ice_frac matrix for the samples using the
% weighting period wt_t_60 (60 days) 
% int = 6; % 6 hr interval to match CCMP winds
% lag = 1:(wt_t *(24/int)); % # observations back in time (1 : # days * # observations per day)
icemat = nan(length(und.time),length(lag)+1); %empty matrix; note: the +1 is to hold the current (i.e. 0 hr lag) observation
t_back = nan(length(und.time),length(lag)+1);
icemat_lat = repmat(und.lat, 1, length(lag)+1);
icemat_lon = repmat(und.lon, 1, length(lag)+1);
X = satice.mdate;
Y = satice.ice_frac;

ice_back = nan(numel(ice_row),numel(lag)+1);
for kk=fliplr(lag)
    t_back(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6hrs ==> lag backwards 720 hrs (30days)
    for m = 1:numel(ice_row)
        ice_back(m,kk) = interpn(satice.mdate,satice.ice_frac(:,ice_row(m),ice_col(m)),t_back(m,kk));
    end
    icemat(:,lag(end)-lag(kk)+1) = ice_back(:,kk); %fill matrix
end

icemat(:,end) = und.ice_frac;
und.icemat = icemat;
und.icemat(isnan(und.icemat)) = 0; % for any nan values within the matrix, set to 0 ice cover;

%--- T, S, MLD and Schmidt number matrices
% Assume constant backwards in time
Tmat = repmat(s_T,1,length(lag)+1);
Smat = repmat(s_S,1,length(lag)+1);
zMLmat = repmat(s_mld,1,length(lag)+1);

[~,Scmat_CH4] = gasmoldiff(Smat,Tmat,'CH4');
[Scmat_N2O] = calc_schmidt(Smat,Tmat,'n2o');

% calculate open water gas transfer velocity following Ho et al. (2006)
spd = 60*60*24; % seconds per day
k.CH4 = kgas(windmat,Scmat_CH4,'Ho06').* spd; % m/d
k.N2O = kgas(windmat,Scmat_N2O,'Ho06') .* spd; % m/d

gd.k_wt_60_CH4 = nan.*s_time;
gd.k_wt_60_N2O = nan.*s_time;

ice = und.icemat;
ice_60 = mean(ice,2); % 60-day mean ice cover

% calculate 60-day weighted gas transfer velocity based on historical ice
% and wind speed data
for kk = 1:length(Scmat_CH4(:,1))
    gd.k_wt_60_CH4(kk) = kw_weighting_ice(k.CH4(kk,:), int/24, wt_t, zMLmat(kk,:),ice(kk,:)./100); % m/d
    gd.k_wt_60_N2O(kk) = kw_weighting_ice(k.N2O(kk,:), int/24, wt_t, zMLmat(kk,:),ice(kk,:)./100); % m/d
end
clear kk k

% calculate 60-day weighted sea-air flux of gases
% Units are umol/m2/d = m/d * umol/kg * kg/m3
gd.F_CH4_60 = gd.k_wt_60_CH4.*(s_ch4-s_ch4_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_60 = gd.k_wt_60_N2O.*(s_n2o-s_n2o_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d

disp('processed 60-day weighted gas fluxes')

gd.slp_inst = slp_inst; % instantaneous slp
gd.slp_30 = slp_30; % 30 day mean slp
gd.ice_inst = ice_inst; % instantaneous ice fraction
gd.ice_30 = ice_30; % 30 day mean ice fraction

% save the data to a new mat file
demo_gas_flux = gd;
save demo_gas_flux.mat demo_gas_flux;


% plot the data to compare concentration, instantaneous, and 60-day fluxes
% axis limits for fluxes
min_F_CH4 = min([gd.F_CH4_inst;gd.F_CH4_60]);
max_F_CH4 = max([gd.F_CH4_inst;gd.F_CH4_60]);

min_F_N2O = min([gd.F_N2O_inst;gd.F_N2O_60]);
max_F_N2O = max([gd.F_N2O_inst;gd.F_N2O_60]);

figure(1)
clf; hold on;
ms=20; % marker size
subplot(2,3,1)
hold on; box on;
scatter(s_lon,s_lat,ms,gd.ch4_nmolkg,'filled','markeredgecolor','k');
xlabel('longitude');
ylabel('latitude');
h=colorbar('southoutside');
xlabel(h, 'CH_4 (nmol kg^{-1})')
title('CH_4')

subplot(2,3,2)
hold on; box on;
scatter(s_lon,s_lat,ms,gd.F_CH4_inst,'filled','markeredgecolor','k');
xlabel('longitude');
ylabel('latitude');
h=colorbar('southoutside');
caxis([min_F_CH4 max_F_CH4]);
xlabel(h, 'CH_4 flux, inst. (\mumol m^{-2} d^{-1})')
title('CH_4 flux, instantaneous')

subplot(2,3,3)
hold on; box on;
scatter(s_lon,s_lat,ms,gd.F_CH4_60,'filled','markeredgecolor','k');
xlabel('longitude');
ylabel('latitude');
h=colorbar('southoutside');
caxis([min_F_CH4 max_F_CH4]);
xlabel(h, 'CH_4 flux, 60-day (\mumol m^{-2} d^{-1})')
title('CH_4 flux, 60-day weighted')

subplot(2,3,4)
hold on; box on;
scatter(s_lon,s_lat,ms,gd.n2o_nmolkg,'filled','markeredgecolor','k');
xlabel('longitude');
ylabel('latitude');
h=colorbar('southoutside');
xlabel(h, 'N_2O (nmol kg^{-1})')
title('N_2O')

subplot(2,3,5)
hold on; box on;
scatter(s_lon,s_lat,ms,gd.F_N2O_inst,'filled','markeredgecolor','k');
xlabel('longitude');
ylabel('latitude');
h=colorbar('southoutside');
caxis([min_F_N2O max_F_N2O]);
xlabel(h, 'N_2O flux, inst. (\mumol m^{-2} d^{-1})')
title('N_2O flux, instantaneous')

subplot(2,3,6)
hold on; box on;
scatter(s_lon,s_lat,ms,gd.F_N2O_60,'filled','markeredgecolor','k');
xlabel('longitude');
ylabel('latitude');
h=colorbar('southoutside');
caxis([min_F_N2O max_F_N2O]);
xlabel(h, 'N_2O flux, 60-day (\mumol m^{-2} d^{-1})')
title('N_2O flux, 60-day weighted')

% display some statistics (optional)
% disp('CH4 flux, instantaneous: mean / std / min / max')
% [mean(gd.F_CH4_inst)  std(gd.F_CH4_inst) min(gd.F_CH4_inst)  max(gd.F_CH4_inst) ]
% 
% disp('CH4 flux, 30 day weighted: mean / std / min / max')
% [mean(gd.F_CH4_30)  std(gd.F_CH4_30) min(gd.F_CH4_30)  max(gd.F_CH4_30) ]
% 
% disp('CH4 flux, 60 day weighted: / std / min / max')
% [mean(gd.F_CH4_60)  std(gd.F_CH4_60) min(gd.F_CH4_60)  max(gd.F_CH4_60) ]
% 
% disp('N2O flux, instantaneous: mean / std / min / max')
% [mean(gd.F_N2O_inst)  std(gd.F_N2O_inst) min(gd.F_N2O_inst)  max(gd.F_N2O_inst) ]
% 
% disp('N2O flux, 30 day weighted: mean / std / min / max')
% [mean(gd.F_N2O_30)  std(gd.F_N2O_30) min(gd.F_N2O_30)  max(gd.F_N2O_30) ]
% 
% disp('N2O flux, 60 day weighted: mean / std / min / max')
% [mean(gd.F_N2O_60)  std(gd.F_N2O_60) min(gd.F_N2O_60)  max(gd.F_N2O_60) ]

     