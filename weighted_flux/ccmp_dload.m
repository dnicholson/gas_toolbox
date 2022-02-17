function ccmp_dload(wdir,fname)

%--------------------------------------------------------------------------
% Download ccmp winds using wget.
% 
% USAGE: ccmp_dload(wdir,fnames);
% 
% INPUT:
% wdir = directory where to save filenames (string)
% fnames = name of text file
% 
% OUTPUT:
% check your wdir directory ;) (there will be files!)
% 
% R. Izett
% UBC Oceanography
% Last modified: Nov. 2016
%--------------------------------------------------------------------------

cd(wdir);

%download the data through command line using wget
    %first, change system directory
        system(['cd ', wdir, '\']); 
    %now download
        system(['wget -i ', fname]); 
    

% % open file with list of nc filenames
%     fid = fopen(fname); 
% 
% %extract the list as a cell array (C)
%     C = textscan(fid,'%s'); C = C{1};
%     
% %close file
%     fclose(fid); clear fid
%     
% %go through each file and download
%     for jj = 1:numel(C);
%         this_file = C{jj};
%         
%         %extract year and month info
%             yr = this_file(20:23);
%             mo = this_file(24:25);
%             
%         %create download url
%             url = strcat('http://data.remss.com/ccmp/v02.0/Y',yr,'/M',mo,'/',this_file);
%         
%         %download through command line using wget
%             %first, change system directory
%                 system(['cd ', wdir, '\']); 
%             %now download
%                 system(['wget -i ', url]); 
%                 
%     end
