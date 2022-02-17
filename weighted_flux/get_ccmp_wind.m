function [mdate, lat, lon, wind] = get_ccmp_wind(wdir,dater,latr,lonr)

%--------------------------------------------------------------------------
% Create list of files for CCMP wind download (.txt file) in specified
% directory.
% 
% USAGE: ccmp_fnames(wdir,sdate,edate)
% 
% INPUT:
% wdir = directory where to save filenames (string)
% dater = start and end dates (matlab datenum format) size [2x1]
% latr = latitude range (-90 to 90) size [2x1]
% lonr = longitude range (-180 to 180) size [2x1]
% 
% OUTPUT:
% mdate = matlab date
% lat = latitude
% lon = longitude
% wind = matrix of wind speed (size: mdate x lat x lon) (units: m/s)
% 
% 
% Based on code by Robert Izett, adapted by Cara Manning
% Last modified: 29 Aug 2021
%--------------------------------------------------------------------------

cd(wdir);
%1) Create file list
    [~] = ccmp_fnames(wdir,dater(1),dater(2)); 
    fname = 'wind_files.txt'; %output from the function is the list of files, but below we just want the name of the .txt file which contains the list

%2) Download files to directory
    wget_dload(wdir,fname); 

%3) Extract wind data from nc file
    [mdate,lat,lon,wind]=ccmp_wind(wdir,fname,latr,lonr); 
    
    clear latr lonr
    
    %plot 1 day
        wi = squeeze(wind(1,:,:));
        figure; pcolor(lon,lat',wi); shading flat
        colorbar;

%4) Delete files from hard drive
    startchar = 44; % number of characters in url preceding the file name
    satfiles_delete(wdir,fname,startchar); % delete the .nc files
    delete([wdir,'\',fname]); % delete the list of filenames
end