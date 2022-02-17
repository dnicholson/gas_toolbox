function [fname] = ccmp_fnames(wdir,sdate,edate)

%--------------------------------------------------------------------------
% Create list of files for CCMP wind download (.txt file) in specified
% directory.
% 
% USAGE: ccmp_fnames(wdir,sdate,edate)
% 
% INPUT:
% wdir = directory where to save filenames (string)
% sdate, edate = start and end dates (matlab date format)
% 
% OUTPUT:
% fname = list of files
% 
% R. Izett
% UBC Oceanography
% Last modified: Nov. 2016
%--------------------------------------------------------------------------

cd(wdir);

%--- specify files you want to download
    %NOTE: the filename format:
    %       CCMP_Wind_Analysis_yyyymmdd_V02.0_L3.0_RSS.nc
    %NOTE ALSO: the download url format:
    %       http://data.remss.com/ccmp/v02.0/Yyyyy/Mmm/CCMP_Wind_Analysis_yyyymmdd_V02.0_L3.0_RSS.nc

    %create vector of dates from start to end date (1 day apart)
        dates = sdate:edate;
      
    %create .txt file to save list of file names
        fid = fopen('wind_files.txt','wt');
        
    %extract dates info and create filenames
        for kk = 1:numel(dates)
           yr = datestr(dates(kk),'yyyy' );
           mo = datestr(dates(kk),'mm' );
           da = datestr(dates(kk),'dd');
            
           fil = strcat(['CCMP_Wind_Analysis_',yr,mo,da,'_V02.0_L3.0_RSS.nc']); %specific filename
           fname(kk,:) = strcat('http://data.remss.com/ccmp/v02.0/Y',yr,'/M',mo,'/',fil); %download url
           
           %wirte files to text file
           fprintf(fid,'%s\n',fname(kk,:));
           
        end
        
    %close file
        fclose(fid);
        
end 