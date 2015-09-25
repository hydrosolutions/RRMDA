function [DEM,R] = getDEM(latlim, lonlim, product)

% DOWNLOAD DEM FROM NASA WORLDWIND SERVER
%
% Index 1: SRTM30 with Bathymetry (900m) merged with global ASTER (30m) 
% Index 2: USGS NED 30m 
% Index 3: ScankortElevationsDenmarkDSM 
% Index 4: ScankortElevationsDenmarkDTM
% Index 5: Aster 30m 
% Index 6: SRTM30 with Bathymetry (900m) merged with global ASTER (30m) and USGS NED (10m)
% Index 7: SRTM30 with Bathymetry (900m) merged SRTM3 V4.1 (90m) and USGS NED (10m)
% Index 8: SRTM30 with Bathymetry (900m) merged SRTM3 V4.1 (90m) and USGS NED (30m)
% Index 9: SRTM3 V4.1 
% Index 10: SRTM30 Plus 
% Index 11: USGS NED 10m



% Find the layers from the NASA WorldWind server.
layers = wmsfind('nasa.network*elev', 'SearchField', 'serverurl');
layers = wmsupdate(layers);

% Select the 'EarthAsterElevations30m' layer containing SRTM30+ data.
source = layers(product);

[DEM, RA] = wmsread(source, 'Latlim', latlim, 'Lonlim', lonlim,'ImageFormat', 'image/bil');
R = refmatToGeoRasterReference(RA,size(DEM));


end
