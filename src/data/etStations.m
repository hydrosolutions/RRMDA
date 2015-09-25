function ETpot_d = etStation(sVal, varargin)

% initial code from Haijing Wang
% adapted and transformed into a function by Lucas Beck
% 19. November 2012

% INPUT VARIABLE sVal needs to be defined as 

% sVal.TEMP -> Temperature in Fahrenheit
% sVal.DEWP -> Dew point (Fahrenheit)
% sVal.WDSP -> Wind speed (knots)
% sVal.MAX -> Maximum Temperature (Fahrenheit)
% sVal.MIN -> Minimum Temperature (Fahrenheit)
% sVal.PRCP -> Precipitation in inch

% OPTIONS for varargin are:
% 'metric' - to have all dimensions in mm (instead of inches), m/s (for
%            knots), and ?C instead of Fahrenheit.
% 'RHUM'   - sVal.RHUM, relative humidity (% [1..100]) instead of sVal.DEWP
%
% 'tmean'  - to do all the calculations based on the mean temperature
%            sVal.TEMP. If not activated then sVal.MAX and sVal.MIN need to 
%            be defined.


if cell2mat(strfind(varargin, 'metric')) == 1
    metric = 1;
else metric = 0;
end

if cell2mat(strfind(varargin, 'RHUM')) == 1
    RHUM = 1;
    try
        sVal.RHUM;
    catch err
        if (strcmp(err.identifier,'MATLAB:nonExistentField'))
            
            msg = sprintf('relative humidity (RHUM) not defined in input');
            error('MATLAB:myCode:missVar', msg);
            % Display any other errors as usual.
        else
            rethrow(err);
        end
        
    end
else RHUM = 0;
end

if cell2mat(strfind(varargin, 'tmean')) == 1
    tm = 1;
else tm = 0;
    try
        sVal.MAX; sVal.MIN;
    catch err
        if (strcmp(err.identifier,'MATLAB:nonExistentField'))
            
            msg = sprintf('MAX and MIN temperature need to be defined');
            error('MATLAB:myCode:missVar', msg);
            
            % Display any other errors as usual.
        else
            rethrow(err);
        end
        
    end
end





%This code is to prepare all the parameters and calculate daily ETpot with
%FAO Penmann equation. n/N is read in from monthly cloudiness data. The
%daily ETpot calculated is writen in ETpot_d.txt.

%keyboard
LON = sVal.LONLAT(:,1);
LAT = sVal.LONLAT(:,2);
z = sVal.LONLAT(:,3);



% keyboard
nN = sVal.nN; % cloudiness

%calculate Patm, gamma, and phi
Patm=101.3*((293-0.0065*z)/293).^5.26;
gamma=0.665*Patm/1000;
phi=pi/180*LAT;

% check if more than one timestep

for fnr=1:size(sVal.LONLAT,1)
    
    for j = 1 : size(sVal.YEARMODA,1)
    % calculate the total number of days within the year (long/short year)
    if isnan(sVal.YEARMODA(j,fnr)) 
        yearDay(j,fnr) = NaN;
        numDay(j,fnr) = NaN;
        Dmonth(j,fnr) = NaN;
    else
    yearDay(j,fnr) = - datenum(strcat( datestr(datenum(num2str(sVal.YEARMODA(j,fnr)), 'yyyymmdd'),'yyyy'), '0101'), 'yyyymmdd') + ...
        datenum(strcat( datestr(datenum(num2str(sVal.YEARMODA(j,fnr)), 'yyyymmdd'),'yyyy'), '1231'), 'yyyymmdd') + 1;
    % the number of days up to date within the current year
    numDay(j,fnr) = - datenum(strcat( datestr(datenum(num2str(sVal.YEARMODA(j,fnr)), 'yyyymmdd'),'yyyy'), '0101'), 'yyyymmdd') ...
        + datenum(num2str(sVal.YEARMODA(j,fnr)), 'yyyymmdd') + 1;
    % calculate also the month for later use
    Dmonth(j,fnr) = str2num(datestr(datenum(num2str(sVal.YEARMODA(j,fnr)),'yyyymmdd'),'mm'));
    end
    end
end

%keyboard
% read in the values passed by sVal:
if metric
    
    for fnr=1:size(sVal.LONLAT,1)
        
        Tmean(:,fnr)= sVal.TEMP(:,fnr); % mean temperature in Celcius
        try
            Tdew(:,fnr)=  sVal.DEWP(:,fnr); % dew point in Celcius
            Tdew(Tdew > 900) = NaN;
        catch
            relhum(:,fnr) = sVal.RHUM(:,fnr); % relative humidity in %
        end
        WDSP(:,fnr)=  sVal.WDSP(:,fnr); %  m/s
        Prcp(:,fnr)=  sVal.PRCP(:,fnr); %  milimeters
        try
            Tmax(:,fnr)=  sVal.MAX(:,fnr); % in Celcius
            Tmin(:,fnr)=  sVal.MIN(:,fnr); % in Celcius
        catch
        end
    end
    
else
    for fnr=1:size(sVal.LONLAT,1)
        
        Tmean(:,fnr)=(sVal.TEMP(:,fnr) - 32)*5/9; % Fahrenheit - Celcius conversion
        try

            Tdew(:,fnr)=(sVal.DEWP(:,fnr) - 32)*5/9; % Fahrenheit - Celcius conversion
            Tdew(Tdew > 900) = NaN;
        catch
            relhum(:,fnr) = sVal.RHUM(:,fnr); % relative humidity in %
        end
        WDSP(:,fnr)= sVal.WDSP(:,fnr) * 0.51444; % knot to m/s
        Prcp(:,fnr)= sVal.PRCP(:,fnr) * 25.4; % from inch to milimeters
        try
            Tmax(:,fnr)=(sVal.MAX(:,fnr) - 32)*5/9; % Fahrenheit - Celcius conversion
            Tmin(:,fnr)=(sVal.MIN(:,fnr) - 32)*5/9; % Fahrenheit - Celcius conversion
        catch
        end
        
        
        
    end
end
    
    
    
    % Tmean=(Tmax+Tmin)/2;
    
    DDelta=4098*0.6108*exp(17.27*Tmean./(Tmean+237.3))./((Tmean+237.3).^2);
    
    if tm == 1
        % calculate es if only Tmean is available:
        es = 0.6108*exp(17.27*Tmean./(Tmean+237.3));
    else
        % calculate es based on Tmax and Tmin
        e0max=0.6108*exp(17.27*Tmax./(Tmax+237.3));
        e0min=0.6108*exp(17.27*Tmin./(Tmin+237.3));
        es=(e0max+e0min)/2;
    end
    
    
    if RHUM == 1
        % ea based on relative humidity
        ea = relhum .* es/100;
        
    else
        % ea based on Dew point
        ea=0.6108*exp(17.27.*Tdew./(Tdew+237.3));
    end
    
    
    dr=1+0.033*cos(2*pi*numDay./yearDay);
    delta=0.409*sin(2*pi*numDay./yearDay-1.39);
%keyboard
for fnr=1:size(sVal.LONLAT,1)

            omegas=acos(-tan(phi(fnr))*tan(delta(:,fnr)));
            Ra=24*60/pi*0.0820*dr(:,fnr).*(omegas*sin(phi(fnr)).*...
                sin(delta(:,fnr))+cos(phi(fnr))*cos(delta(:,fnr)).*sin(omegas));
            % extract the month
            
            for k = 1 : size(Dmonth(:,fnr),1)
                if isnan(Dmonth(k,fnr))
                    Rs(k) = NaN;
                else
                    Rs(k)=(0.25+0.50*nN(fnr,Dmonth(k,fnr)))*Ra(k);
                end
            end
            %keyboard
            Rs0=(0.75+2/10^5*z(fnr))*Ra;
            Rns=(1-0.23)*Rs;
            if tm == 1
                % calculate Rn1 (longwave radiation) based on TEMP (mean
                % temperature)
                Rnl=4.903 / 10^9 * (Tmean(:,fnr)+273).^4 .* ...
                    (0.34-0.14*sqrt(ea(:,fnr))).*(1.35*Rs'./Rs0-0.35);
            else
                Rnl=4.903/10^9*((Tmax(:,fnr)+273).^4+(Tmin(:,fnr)+273).^4)/2 .* ...
                    (0.34-0.14*sqrt(ea(:,fnr))).*(1.35*Rs'./Rs0-0.35);                
            end
            
            Rn=Rns'-Rnl;

            ETpot_d(:,fnr)=(0.408*DDelta(:,fnr).*(Rn-0)+gamma(fnr)*900 .* ...
                WDSP(:,fnr).*(es(:,fnr)-ea(:,fnr))./(Tmean(:,fnr)+273))./...
                (DDelta(:,fnr)+gamma(fnr)*(1+0.34*WDSP(:,fnr)));
    
       
    %end
end



































