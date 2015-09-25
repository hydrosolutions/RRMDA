function cloudiness = cCloud(LON,LAT, cloud)

% calculates the long term average cloudiness of a certain location (LON LAT) and
% returns monthly values
% over open sea, cloudiness is set to 0 (zero) according to original data

% Lucas Beck, 20. Nov 2012

%keyboard


% Dimensions of global cloud matrix are 720 x 360 x 12

% LON-values in matrix are:
lonM = round(359.5 + (2 * LON))+1;

% check and correct limits if necessary
lonM(lonM < 1) = 1;
lonM(lonM > 720) = 720;

% the same for the LAT values
latM = round(179.5 + (2 * LAT))+1;

% check and correct limits if necessary
latM(latM < 1) = 1;
latM(latM > 360) = 360;

% flip up/down cloud data
for i = 1:12
    tCloud(:,:,i) = flipud(cloud(:,:,i)); 
end
cloud = tCloud;

% access global cloud-values
for i = 1 : length(latM)
    cloudiness(i,1:12) = cloud(latM(i), lonM(i), :)/100;
end


























