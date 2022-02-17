function [dat]=extract_ncep(ndir,t0,tf,latr,lonr,type)

%----------------------------------------------------------------------------
%%% ABOUT %%
% This function extracts the time, lat, lon, and gridded data from 
% downloaded ncep/ncar reanalysis 1 files. Files were downloaded using the 
% function download_ncep().
% (https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html)
% 
% [dat]=extract_ncep(ndir,t0,tf,latr,lonr,type);
% 
% INPUT:
%     ndir = directory where downloaded files are saved
%     t0, tf = start / end dates (matlab date format)
%         NOTE: t0 and tf must be within the SAME year
%     latr/lonr = lat/lon ranges ([min max]) 
%         NOTE: longitude = deg E 
%     type = string specifying data you wish to extract; one of: 
%           'wind' (net wind speed),'air' (air temperature), 'slp' (sea
%           level pressure),'rad' (latent heat, sensible heat, net short
%           wave and net long wave radiation), 'rhum' (relative humidity),
%           'prate' (precipitation rate)
%
% OUTPUT:
%     dat = structure containing:
%       mdate = array of matlab times
%       lat, lon = lat/lon arrays
%       units = units of extracted data
%       data = data structure containing specified data type
%           data size: mdate x lat x lon
% 
% R. Izett (rizett{at}eoas.ubc.ca)
% UBC Oceanography
% Last modified: July 2019
% Modified by Cara Manning 29 Aug 2021
%--------------------------------------------------------------------------

%--- CD to specified folder where data are saved 
    %cd([ndir,'\',num2str(year(t0))])
    
%--- get the filename(s) for the specified data type
    if strcmp(type,'wind'); %u and v wind components loaded separately
        fname(1) = dir('uwnd*.nc');
        fname(2) = dir('vwnd*.nc');
        
        %variable names in netcdf files
        var(1,:) = 'uwnd';
        var(2,:) = 'vwnd';
        
        %units
        dat.units.spd = 'm/s';
        dat.units.dir = 'deg';
        
    elseif strcmp(type,'rad') %all rad components loaded separately
        fname(1) = dir('lhtfl*.nc');
        fname(2) = dir('shtfl*.nc');
        fname(3) = dir('nswrs*.nc');
        fname(4) = dir('nlwrs*.nc');
        
        %variable names in netcdf files
        var(1,:) = 'lhtfl';
        var(2,:) = 'shtfl';
        var(3,:) = 'nswrs';
        var(4,:) = 'nlwrs';
        
        %units
        dat.units.lhtfl = 'W/m2, (+) out of planet/ocean'; 
        dat.units.shtfl = 'W/m2, (+) out of planet/ocean';
        dat.units.nswrs = 'W/m2, (+) out of planet/ocean'; 
        dat.units.nlwrs = 'W/m2, (+) out of planet/ocean';
        
    else %everything else
        fname(1) = dir([type,'*.nc']);
        
        %variable names in netcdf files
        var(1,:) = type;
        
        %units
            if strcmp(type,'air')
                dat.units.air = 'C';
            elseif strcmp(type,'slp')
                dat.units.slp = 'mbar';
            elseif strcmp(type,'rhum')
                dat.units.rhum ='%';
            elseif strcmp(type,'prate')
                dat.units.prate = 'kg/m2/s';
            end    
    end
    
%--- load the file, and extract the data
    for kk = 1:numel(fname)
        if kk == 1; %get time, lat, lon only on first iteration
            %time
                t_start = datenum(str2num(datestr(t0,'yyyy')),1,1); %first day of the year
                
                time = double(ncread(fname(kk).name,'time'));
                time = t_start + datenum(0,0,0,(time - time(1)),0,0); %matlab time
                
                ti = find(time >= t0 & time <= tf); %time within t0 and tf range
                dat.mdate = time(ti);
                
            %lat / long
                lat = double(ncread(fname(kk).name,'lat'));
                lon = double(ncread(fname(kk).name,'lon'));
               
                %Adjust lon E/W
                    li = find(lon>180);
                    lon(li) = lon(li)-360;
                    
                %Find lat / lon data within specified latr / lonr
                    la = find(lat >= latr(1) & lat <= latr(2)); 
                    lo = find(lon >= lonr(1) & lon <= lonr(2));

                    dat.lat = lat(la);
                    dat.lon = lon(lo);
                    
                clear li t_start
        end
        
        %Get data within specified time / lat / lon ranges
            dat.(var(kk,:)) = double(ncread(fname(kk).name,var(kk,:),[lo(1), la(1), ti(1)], [length(lo), length(la), length(ti)]));
            dat.(var(kk,:)) = permute(dat.(var(kk,:)),[3,2,1]);
            
        %unit conversions
            if strcmp(type,'slp')
                dat.(var(kk,:)) = dat.(var(kk,:)) .* 0.01; %Pa --> mbar
            elseif strcmp(type,'air')
                dat.(var(kk,:)) = dat.(var(kk,:)) - 273.15; 
            end
    end
    
%--- Calculate net wind speed & direction if type = 'wind'
    if strcmp(type,'wind')
        %--- calculate wind direction
            [th,r] = cart2pol(dat.uwnd,dat.vwnd);
            th = th .* 180 ./ pi;
            th = th + 90;
            th(th<0) = th(th<0) + 360;
            dat.wind_dir = th;
            clear th r
    
        %--- net speed 
            dat.wind_spd=sqrt(dat.vwnd.^2 + dat.uwnd.^ 2);
        
    end
    
disp(' ');
disp('Data extracted');
dat
    
end
   
