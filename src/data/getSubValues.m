function [out] = getSubValues(areas,sub,data,R)
 


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


 %for each subcatchment       
 for d = 1 : l
    
     % first cut the region of interest
     X = [sub(d).X; sub(d).Y]; % polygon
     XMin = floor(min(X(1,:))); YMin = floor(min(X(2,:)));
     XMax = ceil(max(X(1,:))); YMax = ceil(max(X(2,:)));
     
     reg = [[XMin;XMax;XMax;XMin;XMin;NaN],[YMax;YMax;YMin;YMin;YMax;NaN]];
     [inB,onB] = inpolygon(lonC,latC, reg(:,1),reg(:,2));
   
     
     %3D or 4D matrix: calculate subcatchment mean
     if ndims(data)== 4
         for j=1:size(data,4)
             for i=1:size(data,3)
                 dati=data(:,:,i,j);
                 out(i,d,j) = nansum(areas{d}.*dati(inB))/nansum(areas{d});
             end 
         end
    
     
     else
         for i=1:size(data,3)
             dati=data(:,:,i);
             out(i,d) = nansum(areas{d}.*dati(inB))/nansum(areas{d});
         end
     end
 end
 
 out = out;
           