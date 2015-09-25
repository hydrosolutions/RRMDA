function [iE iC] = findNeighbor(thres, R, DEM)



%% IDENTIFICATION
%Identify the cells which need correction

iE=find(DEM>thres);

% Find nearest cell to each of the correction cells
%grid centroid location
lat = zeros(R.RasterSize(1,1),1); %latitude grid coordinates vector
    for i=1:R.RasterSize(1,1)
        lat(i,1) = R.Latlim(1,2) - i*R.CellExtentInLatitude + R.CellExtentInLatitude/2 ;
    end
lon = zeros(R.RasterSize(1,2),1); %longitude grid coordinates vector
    for i=1:R.RasterSize(1,2)
        lon(i,1) = R.Lonlim(1,1) + i*R.CellExtentInLongitude - R.CellExtentInLongitude/2 ;
    end
    
% create coordinate vectors    
[vLon,vLat]=meshgrid(lon,lat);
vLon=vLon(:);
vLat=vLat(:);

% calculate distances, arclen is a sliced vector (requirement for parfor loop)
arclen_sliced=zeros(length(iE),length(vLon));
for l=1:length(iE)
    parfor k=1:length(vLon)
        arclen_sliced(l,k)= distance(vLat(iE(l)),vLon(iE(l)),vLat(k),vLon(k));
    end
end

%convert vector arclen_sliced to normal vector
for l=1:length(iE)
    arclen(1+(l-1)*length(vLon):l*length(vLon))=arclen_sliced(l,:);
end;

B = reshape(arclen,[length(vLon) length(iE)]);
% Exclude cell above thresholds
B(iE,:) = 100;

% Find nearest cell
for l=1:length(iE)
    [Value,iC(l)]=min(B(:,l)); 
end
iC=iC';


end