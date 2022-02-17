function download_ncep(year,ndir,list)

%----------------------------------------------------------------------------
%%% ABOUT %%
% This functions downloads NCEP/NCAR reanalysis 1 or 2 (see below) data during 
% the specified year to a specified or default directory. 
% (https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html)
% 
% USAGE: download_ncep(year,ndir,list);
% 
% BEFORE RUNNING, please ensure that you have downloaded wget.exe and 
% moved it to your C:/Windows/System32 directory. Restart Matlab after 
% doing this for the first time. 
%
% INPUT:
%   year = YYYY
%   ndir = desired directory where files will be downloaded
%   list = optional string array to specify which files to download.
%     e.g. {'uwnd';'vwnd';'slp'}
%   	select from:
%       'uwnd' (10 m u-direction wind speed) << Reanalysis 2
%       'vwnd' (10 m v-direction wind speed) << Reanalysis 2
%       'air' (10 m atm. T) << Reanalysis 2
%       'rhum' (relative humidity) << Reanalysis 1
%       'slp' (sea level pressure) << Reanalysis 1
%       'lhtfl' (latent heat flux at surface) << Reanalysis 1
%       'shtfl' (sensible heat flux at surface) << Reanalysis 1
%       'nswrs' (net short wave radiation at surface) << Reanalysis 1
%       'nlwrs' (net long wave radiation at surface) << Reanalysis 1
%       'prate' (precipitation rate) << Reanalysis 1
% 
% OUTPUT:
%   check your NCEP directory ;) (there will be files!)
%   downloaded files:
%       uwnd (10 m u-direction wind speed)
%       vwnd (10 m v-direction wind speed)
%       air temp (10 m atm. T)
%       rhum (relative humidity)
%       slp (sea level pressure)
%       lhtfl (latent heat flux at surface)
%       shtfl (sensible heat flux at surface)
%       nswrs (net short wave radiation at surface)
%       nlwrs (net long wave radiation at surface)
%       prate (precipitation rate)
%
% R. Izett (rizett{at}eoas.ubc.ca)
% UBC Oceanography
% Modified by Cara Manning 29 Aug 2021
%--------------------------------------------------------------------------

if nargin < 3
    list = [];
end

%--- CD to directory where data will be saved
        cd(ndir);
        
        %--- Create a folder for the given year's data    
            % Check if folder for current year exists
            if isempty(dir(num2str(year)))
                mkdir(num2str(year))
            end
        
        %--- CD to specific year's folder   
            if ispc; cd([cd,'\',num2str(year)]);
            else; cd([cd,'/',num2str(year)]); end
            

%--- Check if any of the data files have already been downloaded 
    dload(1:10) = 0; %1 = download file; 0 = don't download
    dlist = {'uwnd';'vwnd';'air';'rhum';'slp';'lhtfl';'shtfl';'nswrs';'nlwrs';'prate'}; %potential list of files to download
    %download files if there is currently no file downloaded, or if the
    %existing file is from the current year and older than 7 days
    dat = dir('uwnd*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1) && dat.datenum < datenum(year+1,1,1)
                dload(1) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(1) = 1; 
        end
    dat = dir('vwnd*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(2) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(2) = 1;
        end
    dat = dir('air*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(3) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(3) = 1;
        end
    dat = dir('rhum*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(4) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(4) = 1;
        end
    dat = dir('slp*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(5) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(5) = 1;
        end
    dat = dir('lhtfl*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(6) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(6) = 1;
        end
    dat = dir('shtfl*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(7) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(7) = 1;
        end
    dat = dir('nswrs*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(8) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(8) = 1;
        end
    dat = dir('nlwrs*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(9) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(9) = 1;
        end
    dat = dir('prate*.nc');
        if ~isempty(dat) 
            if datenum(date) - dat.datenum > 7 && dat.datenum < datenum(year+1,1,1)
                dload(10) = 1; %download
                if ispc; delete([cd,'\',dat.name]); else; delete([cd,'/',dat.name]); end
            end
        else
            dload(10) = 1;
        end
    clear dat
    
%--- Compare recommended files to be downloaded against user-defined list    
    if ~isempty(list)
        for kk = 1:numel(dload);
            si = find(strcmp(list,dlist{kk}));
            if isempty(si)
                %if there is no match between possible list of files to
                %download and user-desired list, then don't download
                dload(kk) = 0; 
            end
        end
    end


%--- Make a list of data files to download 
    fid = fopen('ncep_files.txt','wt');
    %Wind
    if dload(1) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/uwnd.10m.gauss.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis2/surface_gauss/uwnd.10m.gauss.', num2str(year),'.nc']);
    end
    if dload(2) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/vwnd.10m.gauss.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis2/surface_gauss/vwnd.10m.gauss.', num2str(year),'.nc']);
    end
    %Air T
    if dload(3) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/air.sig995.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis/surface/air.2m.gauss.', num2str(year),'.nc'])
    end
    %RH
    if dload(4) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/rhum.sig995.', num2str(year),'.nc']);
        fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis/surface/rhum.sig995.', num2str(year),'.nc']);
    end
    %SLP
    if dload(5) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/slp.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis/surface/slp.', num2str(year),'.nc']);
    end
    %Latent heat
    if dload(6) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/lhtfl.sfc.gauss.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/lhtfl.sfc.gauss.', num2str(year),'.nc']);
    end
    %Sensible heat
    if dload(7) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/shtfl.sfc.gauss.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/shtfl.sfc.gauss.', num2str(year),'.nc']);
    end
    %Shortwave Rad.
    if dload(8) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/nswrs.sfc.gauss.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/nswrs.sfc.gauss.', num2str(year),'.nc']);
    end
    %Longwave Rad.
    if dload(9) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/nlwrs.sfc.gauss.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/nlwrs.sfc.gauss.', num2str(year),'.nc']);
    end
    %Precipitation rate
    if dload(10) == 1;
        % fprintf(fid,'%s\n',['ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/prate.sfc.gauss.', num2str(year),'.nc']);
	fprintf(fid,'%s\n',['ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis/surface_gauss/prate.sfc.gauss.', num2str(year),'.nc']);
    end
    fclose(fid);

%--- Download the data through command line using wget
    if ~all(dload == 0)
        %first, change system directory
        if ispc; system(['cd ', '\']); else; system(['cd ', '/']); end
        %now download
        system(['wget -i ', 'ncep_files.txt']); 
    end
    
%--- Delete list of files 
    delete('ncep_files.txt')

display(' ')
display('**************************')
display(['NCEP/NCAR ' num2str(year) ' up-to-date!'])
display(['.nc files in ' cd ':'])
dir('*.nc')
display('**************************')
