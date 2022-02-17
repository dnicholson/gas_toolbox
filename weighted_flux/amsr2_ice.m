function [mdate,lat,lon,ice_percent,ice_percent_uncert]=amsr2_ice(wdir,fname)

%--------------------------------------------------------------------------
% Extract the time, lat, long, and ice coverage (%) info from ASMR2 .nc files
% Because the data is not gridded, we will extract all coordinates
%
% USAGE: [mdate,lat,lon,ice_frac,ice_frac_uncert]=asmr2_ice(wdir,fname,latr,lonr);
% 
% INPUT:
% wdir = directory where to save filenames (string)
% fname = filename of txt file containing all of the download urls
% latr/lonr = lat/long ranges ([min max]) NOTE: W longitude = (-)
% 
% OUTPUT:
% mdate = time-stamp
% lat, long = lat/long vectors
% ice_percent = matrix of ice coverage percent, range 0-100% (size: mdate x lat x long) (units: m/s)
% ice_percent_uncert = matrix of ice coverage percent uncertainty, range 0-100% (size: mdate x lat x long) (units: m/s)
% 
% Cara Manning
% University of British Columbia
% March 2019
%--------------------------------------------------------------------------

cd(wdir);

% open file with list of URLs
    fid = fopen(fname); 

%extract the list as a cell array (C)
    C = textscan(fid,'%s'); C = C{1};
    
%close file
    fclose(fid); clear fid

%dummy variable to hold date
    mdate = nan(numel(C),1);
   
%go through each url
% Format is:
% ftp://osisaf.met.no/archive/ice/conc_amsr/yyyy/mm/ice_conc_nh_polstere-100_amsr2_yyyymmddHHHH.nc
% currently the first file always becomes corrupted (unsure why) so only looking at files 2 onward
    for jj = 1:numel(C)
        fn = C{jj}; fn = fn(51:end); %get just the filename section of the string
        %extract date info
            yy = str2num(fn(32:35));
            mm = str2num(fn(36:37));
            dd = str2num(fn(38:39));
            hh = str2num(fn(40:41));     
            
        
        %extract data
        %get lat/long (first iteration only)
            if jj == 1
                lat = double(ncread(fn,'lat')); % dimension 6 (2d matrix)
                lon = double(ncread(fn,'lon')); % dimension 7 (2d matrix)
                
                % initialize variable ice_conc
                % dimensions are same as lat multiplied by number of time
                % points
                ice_percent = nan.*lat; % 2d matrix
                ice_percent = repmat(ice_percent,1,1,numel(C));               
                ice_percent = permute(ice_percent,[3,1,2]);
                ice_percent_uncert = ice_percent;
            end
        %get time info
           mdate(jj)=datenum(yy,mm,dd,hh,0,0); 

        %ice fraction
            ice_percent_i= double(ncread(fn,'ice_conc'))./100; 
            ice_percent_uncert_i = double(ncread(fn,'total_uncertainty'));

%        %get net wind speed from each time observations
%            this_spd=nan(length(tt),length(la),length(lo));
        ice_percent(jj,:,:) = ice_percent_i;
        ice_percent_uncert(jj,:,:) = ice_percent_uncert_i;
        
        disp([num2str(jj) '/' num2str(numel(C)) ' complete.']);
    end