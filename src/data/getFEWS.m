function [P,R_P,timeP] = getFEWS(today,bblat,bblon)

% Download latest FEWS Net Africa Precipipation data from IRI. 
%
% @param today   - serial date of the latest day to download
% @param bblat   - vector with latitude bounding box coordinates (Dimension: 1x2)
% @param bblon   - vector with longitude bounding box coordinates (Dimension: 1x2)
% @return P      - 3D matrix of daily Precipitation in mm (Latitude, Longitude, Time)
% @return R_P    - Georeference structure of Precipitation matrix
% @return timeP  - Time vector (serial dates)
%
% Usage:
%                 [P,R_P,timeP] = getFEWS(today,bblat,bblon)
%
% File:           getFEWS.m
%
% Created:        06/12/2012
%
% Last modified:  15/02/2013
%
% Author:         Sebastian Stoll (hydrosolutions ltd.)           
%
% Purpose:        Download latest FEWS Net Africa Precipipation data from IRI. 
%
% Description:    Download latest FEWS Net Africa Precipipation data from IRI. 
%
% Revisions:      NA
%
% Copyright (C) 2013 Sebastian Stoll
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.


%% Read latitude, longitude adn time

tStart = datenum('01-Jan-1960');
tStart = datenum(tStart);

name=['http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.CPC/.FEWS/.Africa/.DAILY/.ARC2/.daily/T/1/' num2str(today-tStart) '.0/RANGEEDGES/X/' num2str(bblon(1)) '/' num2str(bblon(2)) '/RANGEEDGES/Y/' num2str(bblat(1)) '/' num2str(bblat(2)) '/RANGEEDGES/.est_prcp/dods'];
ncid = netcdf.open(name);

lon = netcdf.getVar(ncid,1); lon=double(lon);
lat = netcdf.getVar(ncid,2); lat = double(lat);
tt = ncread(name,'T'); % complete time vector first point is 1 Jan 1983, in the record.

%% Get precipitation data

start    = [0 0 (length(tt)-20)]; % start location 
count    = [length(lon) length(lat) 20]; % bounding box area

data = netcdf.getVar(ncid,3,start(:),count(:)); % this is the precipitation data
data(data<0) = NaN; data = double(data);

for i=1:size(data,3)
    data(:,:,i)=fliplr(data(:,:,i));
end

data=permute(data,[2 1 3]);

for i=1:size(data,3)
    data(:,:,i)=flipud(data(:,:,i));
end


%% Create georeference structure and time vector
R_P = georasterref('RasterSize',size(data(:,:,1)),'Latlim',...
    [min(lat,[],1) max(lat,[],1)+abs(lat(1)-lat(2))], 'Lonlim',...
        [min(lon,[],1) max(lon,[],1)+abs(lon(1)-lon(2))]);



% Create time vector
t = netcdf.getVar(ncid,0,(length(tt)-20),20); % date vector of last 20 days 
t=t+tStart; t=double(t);
ds=t;

P=data;
timeP=ds';


end

