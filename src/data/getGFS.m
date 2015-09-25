function [fcast] = getGFS(first,last,bblat,bblon)

% Function to get the GFS forecast data from NOMADS (next five days)
%
% @param first   - first serial date of timeseries
% @param last    - serial date of latest day to download
% @param bblat   - vector with latitude bounding box coordinates (Dimension: 1x2)
% @param bblon   - vector with longitude bounding box coordinates (Dimension: 1x2)
% @return fcast  - Forecast structure containing
%                   FT = Forecasted daily temperature in Celcius; 4D Matrix(Latitude, Longitude, Time, Ensemblemember)
%                   FP = Forecasted daily precipitation in mm; 4D Matrix (Latitude, Longitude, Time, Ensemblemember)
%                   FmaxT = Forecasted daily maximum temperature in Celcius; 4D Matrix (Latitude, Longitude, Time, Ensemblemember)
%                   FminT = Forecasted daily minimum temperature in Celcius; 4D Matrix (Latitude, Longitude, Time, Ensemblemember)
%                   R_F = Georeference structure
%                   timeF = Time vector (serial dates)
%
% Usage:
%                 [fcast] = getGFS(first,last,bblat,bblon)
%
% File:           getGFS.m
%
% Created:        07/12/2012
%
% Last modified:  15/02/2013
%
% Author:         Sebastian Stoll (hydrosolutions ltd.)           
%
% Purpose:        Download GFS data from NOMADS. 
%
% Description:    Download GFS data from NOMADS. 
%
% Revisions:      NA
%
% Copyright (C) 2015 hydrosolutions ltd. Zurich, Switzerland
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.


%% Create file names

date = linspace(first,last,last-first+1);
day = datestr(date,'yyyymmdd');
for i=1:size(day,1)
    name{i} = ['http://nomads.ncep.noaa.gov:9090/dods/gens/gens' day(i,:) '/gep_all_00z'];
end
% information = ncinfo(name);


%% Read forecast subset

for k=1:length(name)
    ncid = netcdf.open(name{k});
    % read lat long vectors
    lon = netcdf.getVar(ncid,4); lon=double(lon); %1D
    lat = netcdf.getVar(ncid,3); lat = double(lat); %1D

    % read IDs
    varidT=netcdf.inqVarID(ncid,'tmp2m');
    varidTmax=netcdf.inqVarID(ncid,'tmax2m');
    varidTmin=netcdf.inqVarID(ncid,'tmin2m');
    varidP=netcdf.inqVarID(ncid,'apcpsfc');



    %find subset grid cells
    ilon = find(lon >= bblon(1) & lon <= bblon(2)); %find subset grid cells
    ilat = find(lat >= bblat(1) & lat <= bblat(2));

    start    = [min(ilon)-1 min(ilat)-1 0 0]; % start location (subtract one as netCDF is 0-based, whereas matlab is 1-bases)
    count    = [length(ilon) length(ilat) 21 21]; % bounding box area



    % read Tmean, Tmax, Tmin, Precip
    tF = netcdf.getVar(ncid,varidT,start(:),count(:)); % 4D
    tmaxF = netcdf.getVar(ncid,varidTmax,start(:),count(:)); % 4D
    tminF = netcdf.getVar(ncid,varidTmin,start(:),count(:)); % 4D
    pF = netcdf.getVar(ncid,varidP,start(:),count(:)); % 4D
    lon = netcdf.getVar(ncid,4,start(1),count(1)); % 1D
    lat = netcdf.getVar(ncid,3,start(2),count(2)); % 1D



    %Change dimension order from (lon,lat,time,ens) to (lat,lon,time,ens)

    for i=1:size(tF,3)
        for j=1:size(tF,4)
            tF(:,:,i)=fliplr(tF(:,:,i,j));
            tmaxF(:,:,i)=fliplr(tmaxF(:,:,i,j));
            tminF(:,:,i)=fliplr(tminF(:,:,i,j));
            pF(:,:,i)=fliplr(pF(:,:,i,j));
        end 
    end


    tF=permute(tF,[2 1 3 4]);
    tmaxF=permute(tmaxF,[2 1 3 4]);
    tminF=permute(tminF,[2 1 3 4]);
    pF=permute(pF,[2 1 3 4]);

    % Create georeference
    R_F = georasterref('RasterSize',size(tF(:,:,1,1)),'Latlim',...
        [min(lat,[],1) max(lat,[],1)+abs(lat(1)-lat(2))], 'Lonlim',...
            [min(lon,[],1) max(lon,[],1)+abs(lon(1)-lon(2))]);


    tF=tF-273.15;
    tminF=tminF-273.15;
    tmaxF=tmaxF-273.15;


    %Generate daily timesteps
    FmeanT=ones(size(tF,1),size(tF,2),(size(tF,3)-1)/4,size(tF,4)) ;
    FmaxT=ones(size(tF,1),size(tF,2),(size(tF,3)-1)/4,size(tF,4)) ;
    FminT=ones(size(tF,1),size(tF,2),(size(tF,3)-1)/4,size(tF,4)) ;
    FP=ones(size(tF,1),size(tF,2),(size(tF,3)-1)/4,size(tF,4)) ;


     for i=5:4:size(tF,3)
         FmeanT(:,:,(i-1)/4,:)= mean(tF(:,:,i-3:i,:),3);
         FP(:,:,(i-1)/4,:)= sum(pF(:,:,i-3:i,:),3);
         FmaxT(:,:,(i-1)/4,:)= max(tmaxF(:,:,i-3:i,:),[],3);
         FminT(:,:,(i-1)/4,:)= min(tminF(:,:,i-3:i,:),[],3);
     end



    %% Create time vector
    t6 = date(k)+size(FP,3); % last point
    t0 = date(k); % first point
    ds = linspace(t0,t6-1,size(FP,3))';


    %return
    fcast{k}.FT=FmeanT;
    fcast{k}.FmaxT=FmaxT;
    fcast{k}.FminT=FminT;
    fcast{k}.FP=FP;
    fcast{k}.R_F=R_F;
    fcast{k}.timeF=ds;
end



end

