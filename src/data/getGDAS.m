function [T,maxT,minT,R_T,timeT] = getGDAS(first,last,bblat,bblon)

% Download GDAS data from NOMADS between specified first day (max. today-6) and today-1.
%
% @param first   - first serial date of timeseries
% @param last    - serial date of latest day to download
% @param bblat   - vector with latitude bounding box coordinates (Dimension: 1x2)
% @param bblon   - vector with longitude bounding box coordinates (Dimension: 1x2)
% @return T      - 3D matrix of daily temperature in Celcius (Latitude, Longitude, Time)
% @return maxT   - 3D matrix of maximum daily temperature in Celcius (Latitude, Longitude, Time)
% @return minT   - 3D matrix of minimum daily temperature in Celcius (Latitude, Longitude, Time)
% @return SM200  - 3D matrix of soil moisture in 2m soil column (Latitude, Longitude, Time)
% @return ET0    - 3D matrix of potential evapotranspiration (Latitude, Longitude, Time)
% @return R_T    - Georeference structure
% @return timeT  - Time vector (serial dates)
%
% Usage:
%                 [T,maxT,minT,SM200,ET0,R_T,timeT] = getGDAS(first,last,bblat,bblon)
%
% File:           getGDAS.m
%
% Created:        06/12/2012
%
% Last modified:  23/02/2015 (TS), catching situations where 06h/12h/18/00h
% data is not yet available.
%
% Author:         Sebastian Stoll, Tobias Siegfried (hydrosolutions ltd.)
%
% Purpose:        Download GDAS data from NOMADS.
%
% Description:    Download GDAS data from NOMADS.
%
% Revisions:      NA
%
% Copyright (C) 2015 hydrosolutions ltd. Zurich, Switzerland
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.


%% Create filenames
% name components
date = linspace(first,last,last-first+1);
m =(ones(last-first+1,1).*4)';
m(1,last-first+1)=1;
idx([cumsum([1 m(m>0)])]) = 1;
date = date(cumsum(idx(1:find(idx,1,'last')-1)));
date = datestr(date,'yyyymmdd');
time = linspace(0,18,4);
time =repmat(time,1,last-first+1);

% create filenames
for i=1:size(date,1)
    if time(i)<10
        name{i}=['http://nomads.ncep.noaa.gov:9090/dods/fnl/fnl' date(i,:) '/fnlflx_0' num2str(time(i)) 'z'];
    else
        name{i}=['http://nomads.ncep.noaa.gov:9090/dods/fnl/fnl' date(i,:) '/fnlflx_' num2str(time(i)) 'z'];
    end
end

%% Read temperature subset
ncid = netcdf.open(name{1});
lon = netcdf.getVar(ncid,2); lon=double(lon);
lat = netcdf.getVar(ncid,1); lat = double(lat);

ilon = find(lon > bblon(1) & lon < bblon(2)); %find subset grid cells
ilat = find(lat > bblat(1) & lat < bblat(2));

start    = [min(ilon)-1 min(ilat)-1 0]; % start location (subtract one as netCDF is 0-based, whereas matlab is 1-bases)
count    = [length(ilon) length(ilat) 1]; % bounding box area

lon = netcdf.getVar(ncid,2,start(1),count(1)); % 1D
lat = netcdf.getVar(ncid,1,start(2),count(2)); % 1D

%Ini Grids and ID
temp=zeros(length(lon),length(lat),length(name)); %temperature grid
tmax=zeros(length(lon),length(lat),length(name)); %max temperature grid
tmin=zeros(length(lon),length(lat),length(name)); %min temperature grid
% soil10=zeros(length(lon),length(lat),length(name)); %min temperature grid
% potET=zeros(length(lon),length(lat),length(name)); %min temperature grid

varidT=netcdf.inqVarID(ncid,'tmp2m'); %id of temperature
varidTmax=netcdf.inqVarID(ncid,'tmax2m'); %id of temperature
varidTmin=netcdf.inqVarID(ncid,'tmin2m'); %id of temperature
% varidSoil=netcdf.inqVarID(ncid,'soilm0_200cm'); %id of soil moisture
% varidPot=netcdf.inqVarID(ncid,'pevpravesfc'); %id od potential evaporation rate

%% Create GDAS Temperature grids

for i=1:length(name)
    
    %try % obsolete - see below!
        
        ncid=netcdf.open(name{i}); %id of file
        temp(:,:,i)   = netcdf.getVar(ncid,varidT,start(:),count(:)); % 3D
        tmax(:,:,i)   = netcdf.getVar(ncid,varidTmax,start(:),count(:)); % 3D
        tmin(:,:,i)   = netcdf.getVar(ncid,varidTmin,start(:),count(:)); % 3D
        %         soil200(:,:,i)   = netcdf.getVar(ncid,varidSoil,start(:),count(:)); % 3D
        %         potET(:,:,i)   = netcdf.getVar(ncid,varidPot,start(:),count(:)); % 3D
        
        % NOTE: below is silly since we do not want to fill with NaN simply
        % because a run has not yet been made available online!
        
%     catch
%         
%         disp(['Warning: GDAS data not yet availbable:' name{i}(46:end)]);
%         temp(:,:,i)   = NaN * temp(:,:,1); % 3D
%         tmax(:,:,i)   = NaN * tmax(:,:,1); % 3D
%         tmin(:,:,i)   = NaN * tmin(:,:,1); % 3D
%         %         soil200(:,:,i)   = netcdf.getVar(ncid,varidSoil,start(:),count(:)); % 3D
%         %         potET(:,:,i)   = netcdf.getVar(ncid,varidPot,start(:),count(:)); % 3D
%     
%     end
    
end



%Change dimension order from (lon,lat,time) to (lat,lon,time)

for i=1:size(temp,3)
    temp(:,:,i)=fliplr(temp(:,:,i));
    tmax(:,:,i)=fliplr(tmax(:,:,i));
    tmin(:,:,i)=fliplr(tmin(:,:,i));
    %     soil200(:,:,i)=fliplr(soil200(:,:,i));
    %     potET(:,:,i)=fliplr(potET(:,:,i));
end

temp=permute(temp,[2 1 3]);
tmax=permute(tmax,[2 1 3]);
tmin=permute(tmin,[2 1 3]);
% soil200=permute(soil200,[2 1 3]);
% potET=permute(potET,[2 1 3]);

% from Kelvin to Celsius
temp=temp-273.15;
tmin=tmin-273.15;
tmax=tmax-273.15;

%Conversion from Watt/m2 to mm potential ET
% T = [-20 -10 0 10 20 30 40];
% l = [2.549 2.525 2.501 2.477 2.453 2.430 2.406];
% l_Tavg = interp1(T,l,temp);
% potET=0.0864*potET./l_Tavg;


% Create georeference
R_T = georasterref('RasterSize',size(temp(:,:,1)),'Latlim',...
    [min(lat,[],1) max(lat,[],1)+abs(lat(1)-lat(2))], 'Lonlim',...
    [min(lon,[],1) max(lon,[],1)+abs(lon(1)-lon(2))]);


%Generate daily timesteps
meanT=ones(size(temp,1),size(temp,2),(size(temp,3)-1)/4);
maxT=ones(size(temp,1),size(temp,2),(size(temp,3)-1)/4);
minT=ones(size(temp,1),size(temp,2),(size(temp,3)-1)/4);
% SM200=ones(size(temp,1),size(temp,2),(size(temp,3)-1)/4);
% ET0=ones(size(temp,1),size(temp,2),(size(temp,3)-1)/4);

for i=5:4:size(temp,3)
    meanT(:,:,(i-1)/4)= mean(temp(:,:,i-3:i),3);
    maxT(:,:,(i-1)/4)= max(tmax(:,:,i-3:i),[],3);
    minT(:,:,(i-1)/4)= min(tmin(:,:,i-3:i),[],3);
    %      SM200(:,:,(i-1)/4)= mean(soil200(:,:,i-3:i),3);
    %      ET0(:,:,(i-1)/4)= mean(potET(:,:,i-3:i),3);
end



%% Create time vector
ds = linspace(first,last-1,(size(temp,3)-1)/4)';


%return
T=meanT;
maxT=maxT;
minT=minT;
timeT=ds;
% SM200=SM200;
% ET0=ET0;



end

