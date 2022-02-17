function ccmp_delete(wdir,fname)

%--------------------------------------------------------------------------
% Delete ccmp wind .nc files from hard drive.
% 
% USAGE: ccmp_delete(wdir,fname)
% 
% INPUT:
% wdir = directory where to save filenames (string)
% fname = filename of txt file containing all of the download urls
% 
% OUTPUT:
% check your wdir directory ;) (there won't be files!)
% 
% R. Izett
% UBC Oceanography
% Last modified: Nov. 2016
% Modified by Cara Namning 2021-10-21
%--------------------------------------------------------------------------

cd(wdir)

% open file with list of URLs
    fid = fopen(fname); 

%extract the list as a cell array (C)
    C = textscan(fid,'%s'); C = C{1};
    
%close file
    fclose(fid); clear fid
    
%go through each file and delete
    for bb = 1:numel(C)
        fn = C{bb}; 
        fn = fn(44:end); %get just the filename section of the string
        disp(fn);
        delete(fn)
    end
    
    delete(fname); % also delete the list of files

    