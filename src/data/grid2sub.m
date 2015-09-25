function [out] = grid2sub(sub,data,R1)
% The function to interpolate gridded data on arbitraily shaped - well -
% shapes as defined in shape files.
%
% USAGE: out = downscaleGridPoly(in,dat,lat,lon,attr,R1)
%
% in:  database structure (with indiv. polygons stored in strct arrays)
% dat: gridded data
% lat: vector of latitudes (upper left corner, otherwise adjust code!)
% lon: vetor of longitudes (upper left corner, otherwise adjust code!)
% required).
% attr: field name under which the information is stored
% R1: reference matrix
%
% Tobias Siegfried, 02/06/08
% revision: 07/22/09
% revision 2: 31/03/2011, this is the the umpfed up version for the
% global conflict in international river basins study.
% revision 3: make the parallel option work! 04 / 04 / 2011 => super
% efficient implementation now with 1) poly definition, 2) polybool, 3)
% polyarea - as simple as that!
% revision 4: streamlined some of the code and kickout unessecary overhead.
% Implementing everything without having to loop through indiv. polygons
% for area-weighted calculations. Now many many times more effective.

 
%% General preparation
l = length(sub); % number of shapes

%grid centroid location
lat = zeros(R1.RasterSize(1,1),1); %latitude grid coordinates vector
for i=1:R1.RasterSize(1,1)
    lat(i,1) = R1.Latlim(1,2) - i*R1.CellExtentInLatitude;% + R1.CellExtentInLatitude/2 ;
end
lon = zeros(R1.RasterSize(1,2),1); %longitude grid coordinates vector
for i=1:R1.RasterSize(1,2)
    lon(i,1) = R1.Lonlim(1,1) + i*R1.CellExtentInLongitude;% - R1.CellExtentInLongitude/2 ;
end

[gCX,gCY]=meshgrid(lon,lat);


% do the grid characteristics
E = referenceEllipsoid(7030,'m'); 


tic
for d = 1 : l
    disp(['Calculating d = ' num2str(d) ' out of ' num2str(l)]);
    %% first cut the region of interest
    X = [sub(d).X; sub(d).Y]; % polygon

    
    XMin = floor(min(X(1,:))); YMin = floor(min(X(2,:)));
    XMax = ceil(max(X(1,:))); YMax = ceil(max(X(2,:)));

    
    reg = [[XMin;XMax;XMax;XMin;XMin;NaN],[YMax;YMax;YMin;YMin;YMax;NaN]];
    [inB,onB] = inpolygon(gCX,gCY, reg(:,1),reg(:,2));

    
    gSSX = [gCX(inB) - R1.DeltaLon/2;...
        gCX(inB) - R1.DeltaLon/2;...
        gCX(inB) + R1.DeltaLon/2;...
        gCX(inB) + R1.DeltaLon/2;...
        gCX(inB) - R1.DeltaLon/2;...
        NaN(size(gCX(inB)))];
    gSSY = [gCY(inB) + R1.DeltaLat/2;...
        gCY(inB) - R1.DeltaLat/2;...
        gCY(inB) - R1.DeltaLat/2;...
        gCY(inB) + R1.DeltaLat/2;...
        gCY(inB) + R1.DeltaLat/2;...
        NaN(size(gCY(inB)))];
    gSSX = gSSX(:); gSSY = gSSY(:);
    
    [gSSX, gSSY] = removeExtraNanSeparators(gSSX(:),gSSY(:));
    [gSSX,gSSY] = poly2cw(gSSX,gSSY);
      
    
    % core algo for the intersection areas
    [latInt, lonInt] = polybool('intersection', gSSX, gSSY, X(1,:), X(2,:));
    areasIntersect = areaint(latInt,lonInt,E);

    keyboard
    % get centroids
    [lonInt,latInt] = polysplit(lonInt,latInt);
    intersecCentroidsLatMax = cellfun(@max, latInt); intersecCentroidsLatMin = cellfun(@min, latInt);
    intersecCentroidsLat = mean([intersecCentroidsLatMax intersecCentroidsLatMin],2);
    intersecCentroidsLonMax = cellfun(@max, lonInt); intersecCentroidsLonMin = cellfun(@min, lonInt);
    intersecCentroidsLon = mean([intersecCentroidsLonMax intersecCentroidsLonMin],2);
    % use the inpolygonsCW function for getting the indices. The function
    % assumes that polygons are ordered CW.
    [iGrid,ind] = inpolygonsCW(intersecCentroidsLat,intersecCentroidsLon,gSSX,gSSY);
    ind = [ind{:}]' % ind is now what we need
    
       
    %% result
    
    %2D
    if ndims(data)== 2
        datCut = data(inB');
        out(1,d) = nansum(areasIntersect.*datCut(ind))/nansum(areasIntersect);
    %3D    
    elseif ndims(data)== 3
         for i=1:size(data,3)
             dat1=data(:,:,i);
             datCut = dat1(inB');        
             out(i,d) = nansum(areasIntersect.*datCut(ind))/nansum(areasIntersect);
         end
    %4D
    elseif ndims(data)== 4
        for j=1:size(data,4)
            for i=1:size(data,3)
                dat1=data(:,:,i,j);
                datCut = dat1(inB');        
                out(i,d,j) = nansum(areasIntersect.*datCut(ind));
            end
        end
    end
        
        
       
end

toc

 


 

 

 

 