function satfiles_delete(wdir,fname,startchar)

%--------------------------------------------------------------------------
% Delete satfiles .nc files from hard drive.
% 
% USAGE: satfiles_delete(wdir)
% 
% INPUT:
% wdir = directory where to save filenames (string)
% fname = filename of txt file containing all of the download urls
% fnlength = number of characters in filename
% 
% OUTPUT:
% check your wdir directory ;) (there won't be files!)
% 
% R. Izett
% UBC Oceanography
% Last modified: Nov. 2016
%--------------------------------------------------------------------------

cd(wdir)

% open file with list of ulrs
    fid = fopen(fname); 

%extract the list as a cell array (C)
    C = textscan(fid,'%s'); C = C{1};
    
%close file
    fclose(fid); clear fid
    
%go through each file and delete
    for bb = 1:numel(C)
        fn = C{bb}; fn = fn(startchar:end); %get just the filename section of the string
        delete(fn) 
    end
    