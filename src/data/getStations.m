function sVal = getStations(ll, ur, st, stDate, varargin)
% Download meteorological station data from ftp://ftp.ncdc.noaa.gov/pub/data/gsod/2012/.
%
% @param ll         - lower left corner of extent
% @param ur         - upper right corner of extent
% @param st         - GSOD station list with st.USAF, st.WBAN, st.LAT, st.LON
% @param stDate     - date for the station (end date)
% @param varargin   - numDays: size of the time series. Max is one year (365 or 366)
% @return sVal      - structure with:
%                       ID: Station ID
%                       YEARMODA: date of the data (YYYY,MM,DD)
%                       TEMP: average daily air temperature in degree Fahrenheit
%                       DEWP: daily dew point temperature in degree Fahrenheit
%                       WDSP: average daily wind speed in meters per second
%                       MAX: maximum daily air temperature in degree Fahrenheit
%                       MIN: minimum daily air temperature in degree Fahrenheit
%                       PRCP: daily precipitation sum in millimeters
%                       FRSHTT: storm,fog,snow occurance code
%                       LONLAT: longitude, latitude and elevation of station
%
%
% Usage:
%                 sVal = getStations(ll, ur, st, stDate, varargin)
%
% File:           getStations.m
%
% Created:        20/11/2012
%
% Last modified:  12/03/2013
%
% Author:         Lucas Beck, Sebastian Stoll (hydrosolutions ltd.)
%
% Purpose:        Extract meteorological station data from ftp://ftp.ncdc.noaa.gov/pub/data/gsod/.
%
% Description:    Extract meteorological station data within a certain area specified by the lower left and upper right coordinates (lat, long).
%
% Revisions:      2015/01/06, code rewritten for fixing ftp connection issues, Tobias Siegfried.
%
% Copyright (C) 2013 Lucas Beck, Sebastian Stoll
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.

%%
% URL for ncdc download
ncdcURL = 'ftp://ftp.ncdc.noaa.gov/pub/data/gsod/';

try
    numDays = varargin{1};
catch
    numDays = 1;
end

% convert variables from station list to common unit standard
LAT2 = (st.LAT/1000)';
LON2 = (st.LON/1000)';
ELEV2 = (st.ELEV/10)';

% find stations within lower left (ll) and upper right (ur) corner
latF = find(LAT2 > ll(2) & LAT2 < ur(2));
lonF = find(LON2 > ll(1) & LON2 < ur(1));
station = strcat(st.USAF(intersect(latF, lonF)),'-',st.WBAN(intersect(latF, lonF)));
stn = str2num(cell2mat(st.USAF(intersect(latF, lonF))));
sVal.ID = stn;

% assign latitude, longitude and elevation to "lonlat"
lonlat = [LON2(intersect(latF, lonF)), LAT2(intersect(latF, lonF)),ELEV2(intersect(latF, lonF))] ;

% DOWNLOAD DATA
% Establish ftp connection
ftpO = ftp(ncdcURL(7:23));



for i = 1 : length(station)
    
    dateYears = (str2num(datestr(stDate,'yyyy')) - str2num(datestr(stDate - numDays,'yyyy')));
    nRecFound = 0;
    
    for j = 1 : 1 + dateYears
        
        yy = str2num(datestr(stDate,'yyyy')) - (j-1);
        fName = strcat(station{i},'-',num2str(yy),'.op.gz');
        cd(ftpO,[ncdcURL(24:end) num2str(yy)]);
        
        try
            
            mget(ftpO,fName);
            gunzip(fName);
            temp = importscript(fName(1:end-3));
            delete(fName);
            fileF = 1; nRecFound = nRecFound + 1;
%             disp([fName ' - FILE LOADED!'])
            
        catch
            
%             disp([fName ' - FILE NOT FOUND!'])
            fileF = 0;
            
        end
        
        if j == 1 && fileF
            
            data{i}.STN = temp.STN;
            data{i}.YEARMODA = temp.YEARMODA;
            data{i}.TEMP = temp.TEMP;
            data{i}.DEWP = temp.DEWP;
            data{i}.WDSP = temp.WDSP;
            data{i}.MAX = temp.MAX;
            data{i}.MIN = temp.MIN;
            data{i}.PRCP = temp.PRCP;
            data{i}.FRSHTT = temp.FRSHTT;
            
        elseif j == 1 && ~fileF
            
            data{i}.STN = [];
            data{i}.YEARMODA = [];
            data{i}.TEMP = [];
            data{i}.DEWP = [];
            data{i}.WDSP = [];
            data{i}.MAX = [];
            data{i}.MIN = [];
            data{i}.PRCP = [];
            data{i}.FRSHTT = [];
        
        elseif j > 1 && fileF
            
            data{i}.STN = [data{i}.STN; temp.STN];
            data{i}.YEARMODA = [data{i}.YEARMODA; temp.YEARMODA];
            data{i}.TEMP = [data{i}.TEMP; temp.TEMP];
            data{i}.DEWP = [data{i}.DEWP; temp.DEWP];
            data{i}.WDSP = [data{i}.WDSP; temp.WDSP];
            data{i}.MAX = [data{i}.MAX; temp.MAX];
            data{i}.MIN = [data{i}.MIN; temp.MIN];
            data{i}.PRCP = [data{i}.PRCP; temp.PRCP];
            data{i}.FRSHTT = [data{i}.FRSHTT; temp.FRSHTT];
            
        end
        
    end
    
    % DETECT MISSING OBSERVATIONS

    sVal.YEARMODA(1:numDays,i) = NaN;
    sVal.TEMP(1:numDays,i) = NaN;
    sVal.DEWP(1:numDays,i) = NaN;
    sVal.WDSP(1:numDays,i) = NaN;
    sVal.MAX(1:numDays,i) = NaN;
    sVal.MIN(1:numDays,i) = NaN;
    sVal.PRCP(1:numDays,i) = NaN;
    sVal.FRSHTT(1:numDays,i) = NaN;
    sVal.LONLAT(i,:) = lonlat(i,:);
    
    if nRecFound > 0
        
        [I,ia,ib] = intersect((stDate-numDays+1 : 1 : stDate)',datenum(num2str(data{i}.YEARMODA),'yyyymmdd'));
        sVal.YEARMODA(ia,i) = data{i}.YEARMODA(ib);
        sVal.TEMP(ia,i)  = data{i}.TEMP(ib);
        sVal.DEWP(ia,i)  = data{i}.DEWP(ib);
        sVal.WDSP(ia,i)  = data{i}.WDSP(ib);
        sVal.MAX(ia,i)  = data{i}.MAX(ib);
        sVal.MIN(ia,i)  = data{i}.MIN(ib);
        sVal.PRCP(ia,i)  = data{i}.PRCP(ib);
        sVal.FRSHTT(ia,i)  = data{i}.FRSHTT(ib);
        
    end
    
end

close(ftpO)
% disp('FTP CLOSED - EXITING')
% Done













