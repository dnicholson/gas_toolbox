function [mdate,lat,lon,wind]=ccmp_wind(wdir,fname,latr,lonr)

%--------------------------------------------------------------------------
% Extract the time, lat, lon, and wind speed (m/s) info from ccmp .nc files
% 
% USAGE: [time,lat,lon,wind]=ccmp_wind(wdir,fname,latr,lonr);
% 
% INPUT:
% wdir = directory where to save filenames (string)
% fname = filename of txt file containing all of the download urls
% latr/lonr = lat/lon ranges ([min max]) NOTE: W longitude = (-)
% 
% OUTPUT:
% mdate = time-stamp
% lat, lon = lat/lon vectors
% wind = matrix of wind speed (size: mdate x lat x lon) (units: m/s)
% 
% R. Izett
% UBC Oceanography
% Last modified: Nov. 2016
% Modified by Cara Manning 29 Aug 2021
%--------------------------------------------------------------------------

cd(wdir);

% open file with list of URLs
    fid = fopen(fname); 

%extract the list as a cell array (C)
    C = textscan(fid,'%s'); C = C{1};
    
%close file
    fclose(fid); clear fid

%dumby variable to hold date
    mdate = [];
   
%go through each url
    for jj = 1:numel(C)
%         fn = C{jj}; fn = fn(44:end); %get just the filename section of the string
%         
%         %extract date info
%             yy = str2num(fn(20:23));
%             mm = str2num(fn(24:25));
%             dd = str2num(fn(26:27));

         fn = C{jj}; fn = fn(44:end); %get just the filename section of the string
         
         %extract date info
             yy = str2num(fn(20:23));
             mm = str2num(fn(24:25));
             dd = str2num(fn(26:27));
  
%         fn = C{jj}; fn = fn(43:end); %get just the filename section of the string
%         
%         %extract date info
%             yy = str2num(fn(23:26));
%             mm = str2num(fn(27:28));
%             dd = str2num(fn(29:30));   
              
        %open nc file
    %        ccmp = ncdataset(fn);
        disp('fn')
        disp(fn)
           tt=double(ncread(fn,'time')); % dimension 3
        %extract data
        %get lat/lon (first iteration only)
        
            if jj == 1
                lat = double(ncread(fn,'latitude')); % dimension 2
                lon = double(ncread(fn,'longitude')); % dimension 1
                li = find(lon>180);
                lon(li) = lon(li)-360;
                clear li
                
                la = find(lat >= latr(1) & lat <= latr(2)); %find data within lat range
                lo = find(lon >= lonr(1) & lon <= lonr(2)); %find data within long range
                
                lat = lat(la);
                lon = lon(lo);
                
                % initialize wind variable
                wind = nan(numel(mdate)*numel(tt),numel(la),numel(lo));
            end
        %get time info
           % tt=double(ncread(fn,'time')); % dimension 3
            tt=tt-tt(1); %tt is the hours at which observations are made (e.g. obs. at 0hr, 64h, 12hr, 18hr of the day)
            mdate=[mdate; datenum(yy,mm,dd,tt,0,0)]; %create date vector
            
        %u and v components of wind
            uwnd=double(ncread(fn,'uwnd')); %dimensions lon x lat x time
            uwnd = permute(uwnd,[3,2,1]);  % permute to dimensions time x lon x lat
            uwnd=uwnd(:,la,lo); %keep only data w/in lat/long ranges
            
            vwnd=double(ncread(fn,'vwnd')); 
            vwnd = permute(vwnd,[3,2,1]);  % permute to dimensions time x lon x lat
            vwnd=vwnd(:,la,lo); %keep only data w/in lat/long ranges
        
        %get net wind speed from each time observations
            this_spd=nan(length(tt),length(la),length(lo));
        
            for kk = 1:numel(tt)
               this_u = squeeze(uwnd(kk,:,:));
               this_v = squeeze(vwnd(kk,:,:)); 
               
               this_spd(kk,:,:) = sqrt(this_u.^2 + this_v.^2);
               
            end
           
        wind((jj-1)*numel(tt)+1:(jj)*numel(tt),:,:) = this_spd;     
            
        disp([yy])
        disp(datestr(mdate(end))) 
        disp([num2str(jj) '/' num2str(numel(C)) ' complete.']);
    end

