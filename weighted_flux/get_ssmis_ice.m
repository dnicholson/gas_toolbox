function [mdate, lat, lon, ice_frac] = get_ssmis_ice(wdir,dater,hemi)

%--------------------------------------------------------------------------
% Create list of files for SSMIS ice download (.txt file) in specified
% directory.
% 
% USAGE: get_ssmis_ice(wdir,dater,hemi)
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
% ice = matrix of fractional ice coverage (size: mdate x lat x long) (value 0-1)
% 
% 
% Cara Manning
% Last modified: 29 Aug 2021
%--------------------------------------------------------------------------

cd(wdir);
%1) Create file list
    % save file list to file
    flist = 'ice_files.txt'; %.txt file which contains the list
    ssmis_fnames(wdir,dater(1),dater(2),hemi,flist);
    disp('done step 1 - ssmis_fnames');

%2) Download files to directory
    wget_dload(wdir,flist);
    disp('done step 2 - wget_dload');

%3) Extract ice data from nc file
     [mdate,lat,lon,ice_frac]=ssmis_ice(wdir,flist); 
     disp('done step 3 - ssmis_ice');
     
%4) Delete files from hard drive
    startchar = 46; % character number in URL where filename starts
    satfiles_delete(wdir,flist,startchar)
    delete([wdir,'\',flist]); % delete the file list    
    disp('done step 4 - satfiles_delete');

end