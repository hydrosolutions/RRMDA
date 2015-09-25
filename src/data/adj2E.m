function [cData] = adj2E(thres, data, R, param, option, DEM, R_DEM)



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

% calculate distances
c=0;
for l=1:length(iE)
    for k=1:length(vLon)
        c=c+1;
        [arclen(c),az(c)] = distance(vLat(iE(l)),vLon(iE(l)),vLat(k),vLon(k));
    end
end
B = reshape(arclen,[length(vLon) length(iE)]);
% Exclude cell above thresholds
B(iE,:) = 100;

% Find nearest cell
for l=1:length(iE)
    [Value,iC(l)]=min(B(:,l)); 
end
iC=iC';


%% ADJUST METEOROLOGICAL DATA
%Calculate Correction Factor


switch option
    case 1 %Polynominal factor correction
        cF=1+param{1}.k1*(DEM(iE)./1000-DEM(iC)./1000)+param{1}.k2*(DEM(iE)./1000-DEM(iC)./1000).^2;
        % Adjust cells above thresholds
        for l=1:size(data,3)
            lP=data(:,:,l);
            lP(iE)=lP(iC).*cF;
            cData(:,:,l)=reshape(lP,size(data,1),size(data,2));
        end
           
    case 2 %Additive lapse rate
        cF=(DEM(iE)-DEM(iC)).*param{2}.lapse/100;
        % Adjust cells above thresholds
        for l=1:size(data,3)
            lT=data(:,:,l);
            lT(iE)=lT(iC)+cF;
            cData(:,:,l)=reshape(lT,size(data,1),size(data,2));
        end
end
  
        



end