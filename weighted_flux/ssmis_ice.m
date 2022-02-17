function [mdate,lat,lon,ice_frac]=ssmis_ice(wdir,fname)

%--------------------------------------------------------------------------
% Extract the time, lat, long, and ice fraction (0-1) info from SSMIS .nc files
% 
% USAGE: [t,lat,long,ice_frac]=ssmis_ice(wdir,fname,latr,lonr);
% 
% INPUT:
% wdir = directory where to save filenames (string)
% fname = filename of txt file containing all of the download urls
% latr/lonr = lat/long ranges ([min max]) NOTE: W longitude = (-)
% 
% OUTPUT:
% mdate = time-stamp
% lat, long = lat/long vectors
% ice_frac = matrix of ice coverage fraction (range 0-1) (size: mdate x lat x long) 
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
% ftp://osisaf.met.no/archive/ice/conc/yyyy/mm/ice_conc_nh_polstere-100_multi_yyyymmddHHHH.nc
% ftp://osisaf.met.no/archive/ice/conc/2015/05/ice_conc_nh_polstere-100_multi_201506101200.nc
for jj = 1:numel(C)
    fn = C{jj}; fn = fn(46:end); %get just the filename section of the string
    %extract date info
    yy = str2num(fn(32:35));
    mm = str2num(fn(36:37));
    dd = str2num(fn(38:39));
    hh = str2num(fn(40:41));
    
    
    %extract data
    %get lat/lon (first iteration only)
    % initialize variable ice_frac
    % dimensions are same as lat multiplied by number of time
    % points    
    if jj == 1
        lat = double(ncread(fn,'lat')); % dimension 6 (2d matrix)
        lon = double(ncread(fn,'lon')); % dimension 7 (2d matrix)
        ice_frac = nan.*lat; % 2d matrix
        ice_frac = repmat(ice_frac,1,1,numel(C));
        ice_frac = permute(ice_frac,[3,1,2]);        
    end

    %get time info
    mdate(jj)=datenum(yy,mm,dd,hh,0,0);
    
   % import ice fraction IF the nc file has data, otherwise set to NaN
    jj_to_remove = [];
    try
      ice_frac_i= double(ncread(fn,'ice_conc'))./100;   
    catch 
      ice_frac_i = nan.*lat;
      jj_to_remove = [jj_to_remove; jj];
      disp('catch')
      disp(datestr(mdate(jj)));
    end

    % code to add in total_uncertainty; only available on some days
%     % check if 'total_uncertainty' is available
%     ncid = netcdf.open(fn,'nowrite');
% 
%     try
%         VarID = netcdf.inqVarID(ncid,'total_uncertainty');
%         ice_conc_uncert_i = ncread(fn,'total_uncertainty');
%     catch exception
%         if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
%             ice_conc_uncert_i = nan.*ice_frac_i;
%         end
%     end
%     
%     netcdf.close(fn);
%     ice_frac_uncert_i = ice_conc_uncert_i./100;

            
        %get ice fraction
        ice_frac(jj,:,:) = ice_frac_i;
%       ice_frac_uncert(jj,:,:) = ice_frac_uncert_i;
        
        disp([num2str(jj) '/' num2str(numel(C)) ' complete.']);
end
% remove any NC files that did not contain any data
if ~isempty(jj_to_remove)
  disp('jj_to_remove');
  disp(jj_to_remove);
  ice_frac(jj_to_remove,:,:) = [];
  mdate(jj_to_remove,:,:) = [];
end
end
