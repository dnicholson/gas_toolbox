%----------------------------------------------------------------------------
%%% DESCRIPTION %%%
% This script shows an example of how to download and extract NCEP/NCAR
% reanalysis 1 data (SLP), SSMIS and AMSR-2 sea ice data, and CCMP wind
% speed data
%
% You will need to add the GitHub repository gas_toolbox to your matlab path
%
% BEFORE RUNNING, please ensure that you have downloaded wget.exe and 
% moved it to your C:/Windows/System32 directory. Restart Matlab after 
% doing this for the first time. 
%
% For the demonstration script, we will calculate the flux from stations
% on the DBO 2018 cruise in the Bering and Chukchi Seas. 
% latitude range is 62 to 72 
% longitude range -176 to -160 
% samples were collected between July 16, 2018 to July 23, 2018
%----------------------------------------------------------------------------

%--- (OPTIONAL) Clear current MATLAB session
	clear all; close all; clc;
   
%----------------------------------------------------------------------------
%%% INPUT PARAMETERS %%%
%--- Indicate file path where downloaded files will be saved temporarily
        % change to reflect your file structure
    	ndir = 'C:\Users\cmanning\Dropbox\gas_toolbox\weighted_flux\data';

%--- Set desired start/end dates and lat/lon ranges for data
	% NOTE: for data sets spanning multiple years, repeat this sequence for 
    % each year of data (i.e. t0 and tf should be in the same year)
    % Because the weighting scheme uses prior data (typically 30 or 60
    % days)
    % t0 needs to be at least 31 or 61 days prior to the first observation
    % tf must be after the last observation
    % The NCEP SLP has a spacing of 2.5 x 2.5 degrees (CCMP and ice are
    % higher resolution), so find the latitude and longitude ranges of your
    % observations and extend by 2.5 degrees.
    t0 = datenum(2018,5,10); % start date
    tf = datenum(2018,7,24,18,0,0); % end date
    latr = [59 75]; % deg N
    lonr = [-178 -158]; % deg E
    hemi = 'nh'; % hemisphere for sea ice data
 
%%   
    
%--------------------------------------------------------------------------
%%% RUN DOWNLOADING / EXTRACTING SCRIPTS %%%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% download NCEP SLP data
% -------------------------------------------------------------------------

%--- Download files to directory 
    % This function is able to download many different variables. In this
    % case we only need the SLP. Consult the function readme and modify
    % accordingly to add or remove other variables. 
    % The function will not download files that already exist on your hard
    % drive, unless they are > 7 days old and from the current year.
	% Watch as files start to populate your ndir. Matlab will also show the
	% download progress.
    
    % download NCEP SLP data
    download_ncep(str2num(datestr(t0,'yyyy')),ndir,'slp'); 
    
%--- Extract NCEP SLP data
    dat=extract_ncep(ndir,t0,tf,latr,lonr,'slp');

	%plot average to check it worked
        wi = squeeze(nanmean(dat.slp,1));
        figure; 
        pcolor(dat.lon,dat.lat',wi); shading flat;
        xlabel('longitude'); ylabel('latitude'); 
        title([datestr(t0,'yyyy'), ' average slp (mbar)']);       
        colorbar;
    
demo_SLP = dat;

save demo_SLP.mat demo_SLP;

clear dat wi; % clear data to free up memory
%%
% -------------------------------------------------------------------------
% download CCMP wind speed data
% -------------------------------------------------------------------------
dater = [t0 tf];
[mdate, lat, lon, wind] = get_ccmp_wind(ndir,dater,latr,lonr);

demo_CCMP.mdate = mdate;
demo_CCMP.lat = lat;
demo_CCMP.lon = lon;
demo_CCMP.wind = wind;

	%plot average to check it worked
        wi = squeeze(nanmean(demo_CCMP.wind,1));
        figure; 
        pcolor(demo_CCMP.lon,demo_CCMP.lat',wi); shading flat;
        xlabel('longitude'); ylabel('latitude'); 
        title([datestr(t0,'yyyy'), ' average wind speed (m/s)']);       
        colorbar;


save demo_CCMP.mat demo_CCMP;

clear demo_CCMP mdate lat lon wind; % clear data to free up memory

%%
% -------------------------------------------------------------------------
% download SSMIS sea ice data 
% SSMIS data is available from 17 Aug 2009 to present
% -------------------------------------------------------------------------
dater = [t0 tf];
[mdate, lat, lon, ice_frac] = get_ssmis_ice(ndir,dater,hemi);

demo_SSMIS.mdate = mdate;
demo_SSMIS.lat = lat;
demo_SSMIS.lon = lon;
demo_SSMIS.ice_frac = ice_frac;

save demo_SSMIS.mat demo_SSMIS;

clear demo_SSMIS mdate lat lon ice_frac; % clear data to free up memory

%%
% -------------------------------------------------------------------------
% download AMSR-2 sea ice data 
% AMSR-2 data is available from 19 Sep 2016 to present
% -------------------------------------------------------------------------
dater = [t0 tf];
[mdate, lat, lon, ice_frac] = get_amsr2_ice(ndir,dater,hemi);

demo_AMSR2.mdate = mdate;
demo_AMSR2.lat = lat;
demo_AMSR2.lon = lon;
demo_AMSR2.ice_frac = ice_frac;

save demo_AMSR2.mat demo_AMSR2;

clear demo_AMSR2 mdate lat lon ice_frac; % clear data to free up memory
