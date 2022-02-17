function [mdate, lat, lon, ice_frac,ice_frac_uncert] = get_amsr2_ice(wdir,dater,hemi)

%--------------------------------------------------------------------------
% Create list of files for AMSR2 ice download (.txt file) in specified
% directory.
% 
% USAGE: ccmp_fnames(wdir,sdate,edate)
% 
% INPUT:
% wdir = directory where to save filenames (string)
% dater = start and end dates (matlab datenum format), size [2x1]
% hemi = hemisphere of interest ('nh' or 'sh')
% 
% OUTPUT:
% mdate = matlab date
% lat = latitude
% lon = longitude
% ice = matrix of ice coverage fraction (size: mdate x lat x lon) (value 0-1)
% 
% 
% Cara Manning
% UBC 
% Last modified: Mar 2019
%--------------------------------------------------------------------------

cd(wdir);
%1) Create file list
    % save file list to file
    flist = amsr2_fnames(wdir,dater(1),dater(2),hemi); 
    fname = 'ice_files.txt'; %output from the function is the list of files, but below we just want the name of the .txt file which contains the list

%2) Download files to directory
    wget_dload(wdir,fname); 
    
%3) Extract ice data from nc file
     [mdate,lat,lon,ice_frac,ice_frac_uncert]=amsr2_ice(wdir,fname); 

%4) Delete files from hard drive
    fnlength = 51; % character number in URL where filename starts
    satfiles_delete(wdir,fname,fnlength)
    delete([wdir,'\',fname]); % delete the file list
end