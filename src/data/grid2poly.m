function [out] = grid2poly(sub,data,R)
% Function to downscale gridded data (3D or 4D) on arbitraily shaped polygons and returns an area weighted average.
%
% @param sub     - structure with indiv. polygons stored in strct arrays
% @param data    - gridded data (lat,long,time,ensemble)
% @param R       - georeference structure of data
% @return out    - data for each timestep and subcatchment and optionally
%                  ensemble member (Dimesions: time,subcatchment,ensemble member)
% 
% Usage:
%                 [out] = grid2poly(sub,data,R)
%
% File:           grid2poly.m
%
% Created:        19/12/2012
%
% Last modified:  04/12/2014
%
% Author:         Tobias Siegfried, Sebastian Stoll (hydrosolutions ltd.)           
%
% Purpose:        Downscale gridded data on arbitraily shaped polygons and returns an area weighted average.
%
% Description:    Downscale gridded data on arbitraily shaped polygons and returns an area weighted average.
%
% Revisions:      NA
%
% Copyright (C) 2013 Tobias Siegfried, Sebastian Stoll
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.


%% General preparation
l = length(sub); % number of shapes

%grid centroid location
lat = zeros(R.RasterSize(1,1),1); %latitude grid coordinates vector
for i=1:R.RasterSize(1,1)
    lat(i,1) = R.Latlim(1,2) - i*R.CellExtentInLatitude + R.CellExtentInLatitude/2 ;
end
lon = zeros(R.RasterSize(1,2),1); %longitude grid coordinates vector
for i=1:R.RasterSize(1,2)
    lon(i,1) = R.Lonlim(1,1) + i*R.CellExtentInLongitude - R.CellExtentInLongitude/2 ;
end
[lon,lat]=meshgrid(lon,lat);

latC = lat(:)'; % grid centroid location X
lonC = lon(:)'; % grid centroid location Y

earthradius = earthRadius('m');

 %for each subcatchment       
 for d = 1 : l
    
     % first cut the region of interest
     X = [sub(d).X; sub(d).Y]; % polygon
     XMin = floor(min(X(1,:))); YMin = floor(min(X(2,:)));
     XMax = ceil(max(X(1,:))); YMax = ceil(max(X(2,:)));
     
     reg = [[XMin;XMax;XMax;XMin;XMin;NaN],[YMax;YMax;YMin;YMin;YMax;NaN]];
     [inB,onB] = inpolygon(lonC,latC, reg(:,1),reg(:,2));
     
     
     gSSX = [lonC(inB) - abs(lon(1,1)-lon(1,2))/2;...
            lonC(inB) + abs(lon(1,1)-lon(1,2))/2;...
            lonC(inB) + abs(lon(1,1)-lon(1,2))/2;...
            lonC(inB) - abs(lon(1,1)-lon(1,2))/2;...
            lonC(inB) - abs(lon(1,1)-lon(1,2))/2;...
            NaN(size(lonC(inB)))];
            
     gSSY = [latC(inB) + abs(lat(1,1)-lat(2,1))/2;...
            latC(inB) + abs(lat(1,1)-lat(2,1))/2;...
            latC(inB) - abs(lat(1,1)-lat(2,1))/2;...
            latC(inB) - abs(lat(1,1)-lat(2,1))/2;...
            latC(inB) + abs(lat(1,1)-lat(2,1))/2;...
            NaN(size(latC(inB)))];
     
     [gSSX, gSSY] = removeExtraNanSeparators(gSSX(:),gSSY(:));
     [gSSX,gSSY] = poly2cw(gSSX,gSSY); 
     
     [latc,lonc] = polysplit(gSSY,gSSX);
            
     areas = zeros(length(latc),1);
     
     % Calculate corresponding areas
     for i = 1 : length(latc)
         
            [lonInt, latInt] = polybool('intersection', lonc{i}, latc{i}, X(1,:), X(2,:));
            
            if ~isempty(lonInt)
                areas(i,1) = nansum(areaint(latInt,lonInt,earthradius));
            else
                areas(i,1) = 0;
            end
     end
     
     %3D or 4D matrix
     if ndims(data)== 4
         
         for j=1:size(data,4)
             
             for i=1:size(data,3)
                 dati=data(:,:,i,j);
                 out(i,d,j) = nansum(areas.*dati(inB)')/nansum(areas);

             end 
         end
     
     
     else
             
         for i=1:size(data,3)
             dati=data(:,:,i);
             out(i,d) = nansum(areas.*dati(inB)')/nansum(areas);
             
         end
     end
   
 end
 
 out = out;
     
    
     
           