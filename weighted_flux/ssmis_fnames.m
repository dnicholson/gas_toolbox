function ssmis_fnames(wdir,sdate,edate,hemi,flist)

%--------------------------------------------------------------------------
% Create list of files for SSMIS ice download (.txt file) in specified
% directory.
% 
% USAGE: ssmis_fnames(wdir,sdate,edate,hemi,fname)
% 
% INPUT:
% wdir = directory where to save filenames (string)
% sdate, edate = start and end dates (matlab date format)
% hemi = hemisphere; 'N' for northern hemisphere or 'S' for southern hemisphere
% flist = file where list of .nc file URLs is saved
% 
% OUTPUT:
% fname = list of files
% 
% Cara Manning
% Last modified: Feb 2019
%--------------------------------------------------------------------------

cd(wdir);

%--- specify files you want to download
    % NOTE: the filename format:
    %       ice_conc_nh_polstere-100_amsr2_201707011200.nc
    % NOTE ALSO: the download url format:
    % ftp://osisaf.met.no/archive/ice/conc_amsr/yyyy/mm/ice_conc_nh_polstere-100_amsr2_yyyymmddHHHH.nc
    % create vector of dates from start to end date (1 day apart)
        dates = sdate:edate;
    
    % create .txt file to save list of file names
        fid = fopen(flist,'wt');
        
        flist = char.empty(numel(dates),0);
    % extract dates info and create filenames
        for kk = 1:numel(dates)
           yyyy = num2str(year(dates(kk)));
           mm = month(dates(kk));
           dd = day(dates(kk));     
           
           fil = ['ice_conc_',hemi,'_polstere-100_multi_',num2str(yyyy),num2str(mm,'%02d'),num2str(dd,'%02d'),'1200.nc'];  
          flist{kk} = ['ftp://osisaf.met.no/archive/ice/conc/',num2str(yyyy),'/',num2str(mm,'%02d'),'/',fil];
           
           fprintf(fid,'%s\n',flist{kk});
           
        end
        
    % close file
        fclose(fid);
        
end 