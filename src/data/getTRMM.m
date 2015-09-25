function [P,R_P,timeP] = getTRMM(today,bblat,bblon)

% Download latest TRMM precipitation data from IRI. 
%
% @param today   - serial date of the latest day to download
% @param bblat   - vector with latitude bounding box coordinates (Dimension: 1x2)
% @param bblon   - vector with longitude bounding box coordinates (Dimension: 1x2)
% @return P      - 3D matrix of daily Precipitation in mm (Latitude, Longitude, Time)
% @return R_P    - Georeference structure of Precipitation matrix
% @return timeP  - Time vector (serial dates)
%
% Usage:
%                 [P,R_P,timeP] = getTRMM(today,bblat,bblon)
%
% File:           getTRMM.m
%
% Created:        22/01/2015
%
% Last modified:  22/01/2015
%
% Author:         Sebastian Stoll (hydrosolutions ltd.)           
%
% Purpose:        Download latest TRMM Precipipation data from IRI. 
%
% Description:    Download latest TRMM Precipipation data from IRI. 
%
% Revisions:      NA
%
% Copyright (C) 2015 hydrosolutions ltd.
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.


%% Read latitude, longitude adn time

tStart = datenum('01-Jan-1998');
tStart = datenum(tStart);

% number of days to download
nD=100;

name=['http://iridl.ldeo.columbia.edu/SOURCES/.NASA/.GES-DAAC/.TRMM_L3/.TRMM_3B42RT/.v7/.daily/.precipitation/T/1/' num2str(today-tStart) '.0/RANGEEDGES/X/' num2str(bblon(1)) '/' num2str(bblon(2)) '/RANGEEDGES/Y/' num2str(bblat(1)) '/' num2str(bblat(2)) '/RANGEEDGES/dods'];
ncid = netcdf.open(name);

lon = netcdf.getVar(ncid,1); lon=double(lon);
lat = netcdf.getVar(ncid,2); lat = double(lat);
tt = ncread(name,'T'); % complete time vector first point is 1 Jan 1998, in the record.

%% Get precipitation data

start    = [0 0 (length(tt)-nD)]; % start location 
count    = [length(lon) length(lat) nD]; % bounding box area

data = netcdf.getVar(ncid,3,start(:),count(:)); % this is the precipitation data
data(data<0) = NaN; data = double(data);

for i=1:size(data,3)
    data(:,:,i)=fliplr(data(:,:,i));
end

data=permute(data,[2 1 3]);


%% Create georeference structure and time vector
R_P = georasterref('RasterSize',size(data(:,:,1)),'Latlim',...
    [min(lat,[],1) max(lat,[],1)+abs(lat(1)-lat(2))], 'Lonlim',...
        [min(lon,[],1) max(lon,[],1)+abs(lon(1)-lon(2))]);
R_P.ColumnsStartFrom = 'north';  



% Create time vector
t = netcdf.getVar(ncid,0,(length(tt)-nD),nD); % date vector 
t=t+tStart; t=double(t);
ds=t;

P=data;
timeP=ds';


end

