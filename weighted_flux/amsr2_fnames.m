function [fname] = amsr2_fnames(wdir,sdate,edate,hemi)

%--------------------------------------------------------------------------
% Create list of AMSR2 sea ice files for download (.txt file) in specified
% directory.
% 
% USAGE: amsr2_fnames(wdir,sdate,edate)
% 
% INPUT:
% wdir = directory where to save filenames (string)
% sdate, edate = start and end dates (matlab date format)
% hemi = hemisphere; 'N' for northern hemisphere or 'S' for southern hemisphere
% 
% OUTPUT:
% fname = list of files
% 
% Cara Manning 
% Modified: 29 Aug 2021
%--------------------------------------------------------------------------

cd(wdir);

%--- specify files you want to download
    %NOTE: the filename format:
    %       ice_conc_nh_polstere-100_amsr2_201707011200.nc
    
    %NOTE ALSO: the download url format:
    %ftp://osisaf.met.no/archive/ice/conc_amsr/yyyy/mm/ice_conc_nh_polstere-100_amsr2_yyyymmddHHHH.nc
    %create vector of dates from start to end date (1 day apart)
        dates = sdate:edate;
    
    %create .txt file to save list of file names
        fid = fopen('ice_files.txt','wt');
        
        fname = char.empty(numel(dates),0);
    %extract dates info and create filenames
        for kk = 1:numel(dates)
           yyyy = num2str(year(dates(kk)));
           mm = month(dates(kk));
           dd = day(dates(kk));

           % add leading zeros to months and days where needed
 %          if mm < 10
 %              mm = strcat(sprintf('%02d',mm));
 %          end
 %          if dd < 10
 %              dd = strcat(sprintf('%02d',dd));
 %          end        
           
           fil = ['ice_conc_',hemi,'_polstere-100_amsr2_',num2str(yyyy),num2str(mm,'%02d'),num2str(dd,'%02d'),'1200.nc'];  
          fname{kk} = ['ftp://osisaf.met.no/archive/ice/conc_amsr/',num2str(yyyy),'/',num2str(mm,'%02d'),'/',fil];

           %write files to text file
           fprintf(fid,'%s\n',fname{kk});
           
        end
        
    %close file
        fclose(fid);
        
end 