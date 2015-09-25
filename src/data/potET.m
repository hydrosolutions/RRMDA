function [ETo] = potET(tmean,tmax,tmin,lat,time)

%% last modified 23/11/2012 Sebastian Stoll
%
% Function to calculate potential evapotranspiration according to
% Hargreaves
%
% Input:
% tmean = daily mean temperature (Celcius)
% tmax = daily maximum temperature (Celcius)
% tmin = daily minimum temperature (Celcius)
% lat = latitude matrix (2D)
% time = date vector (1D serial dates)
%
% Output:
% ETo = potential evapotranspiration (mm/d)





%% Calculate Potential evapotranspiration according to Hargreaves [mm/day]
% Day of year: 
year = datevec(time);
year(:,2:end) = 0; % set 0 except for year
yearNum = datenum(year);
dy = time - yearNum; %calculate day of the year

% Calculate SolarDeclination vector
delta = asin(0.4*sin((2*pi/365)*(dy-82)));    
    
% Calculate Eccentricity correction of the earth's orbit vector
E0 = 1+0.033*cos(2*pi*dy/365);

% Calculate Latitude (rad) matrix
phi = lat*2*pi/360;

% Calculate angular velocity of the earth's rotation [rad/hr]
omega = 0.2618;


%Conversion to ET (mm) equivalent
T = [-20 -10 0 10 20 30 40];
l = [2.549 2.525 2.501 2.477 2.453 2.430 2.406];
l_Tavg = interp1(T,l,tmean);


% Calculate Potential evapotranspiration 
for i=1:length(delta)
    % Calculate hour of sunrise vector
    hsr(:,:,i)= acos(-tan(delta(i)).*tan(phi))./omega;  
    % Calculate extraterrestrial radiation
    H0(:,:,i)=37.59*E0(i).*(omega*hsr(:,:,i).*(sin(delta(i)).*sin(phi))+(cos(delta(i)).*cos(phi).*sin(omega*hsr(:,:,i))));
   
    %3D or 4D matrix
    if ndims(tmean)== 4
        for j=1:size(tmean,4)
            lETo(:,:,i,j)= (0.0023.*H0(:,:,i).*(tmax(:,:,i,j)-tmin(:,:,i,j)).^0.5).*(tmean(:,:,i,j)+17.8);
        end
    else
        lETo(:,:,i)= (0.0023.*H0(:,:,i).*(tmax(:,:,i)-tmin(:,:,i)).^0.5).*(tmean(:,:,i)+17.8); 
    end
end  
ETo=lETo./l_Tavg;


