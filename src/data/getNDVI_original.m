function [A,set,R] = getNDVI(inDate,region)
% Download latest eNDVI image from http://earlywarning.usgs.gov/ftp2/africa/emodis/. 
%
% @param inDate     - wished date
% @return A         - eNDVI image
% @return set       - number of image in current year
% 
% Usage:
%                 [A,set] = getNDVI(inDate,region)
%
% File:           getNDVI.m
%
% Created:        30/11/2012
%
% Last modified:  01/12/2014
%
% Author:         Lucas Beck, Sebastian Stoll (hydrosolutions ltd.)           
%
% Purpose:        Download latest eNDVI image from http://earlywarning.usgs.gov/ftp2/africa/emodis/. 
%
% Description:    Download latest eNDVI image from http://earlywarning.usgs.gov/ftp2/africa/emodis/. 
%
% Revisions:      NA
%
% Copyright (C) 2013 Lucas Beck, Sebastian Stoll
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.


%% eNDVI download


switch region
    case 'North Africa'
        reg='africa';
        regFile='na';
        
    case 'West Africa'
        reg='africa';
        regFile='wa';
        
    case 'East Africa'
        reg='africa';
        regFile='ea';
        
    case 'South Africa'
        reg='africa';
        regFile='sa';
        
    case 'Central America'
        reg='centralamerica';
        regFile='ca';
        
    case 'Central Asia'
        reg='centralasia';
        regFile='cta';
        
    case 'Yemen'
        reg='yemen';
        regFile='yem';        
end


% location of NDVI file
fullNDVI = strcat('http://earlywarning.usgs.gov/ftp2/',reg,'/emodis/');

% create number for actual period -> set
set = ceil((inDate - datenum((strcat('0101',datestr(inDate,'yyyy'))),'ddmmyyyy'))/5);

sFlag = 0; % flag for while loop
while sFlag == 0
    if set < 0
        setStr = (['0',num2str(set)]);
    else
        setStr = num2str(set);
    end 
        filename = strcat(regFile, setStr, datestr(inDate,'yy'),'.zip');
        % try to access the file. If not exist, try to access file from
        % earlier set.
        try
            unzip(filename);
%             disp([filename, ' already exists. Latest data up to date']);
            sFlag = 1;
            
        catch
%             disp(['trying to download ', filename]);
            [filestr,status] = urlwrite(strcat(fullNDVI,filename),filename);
            if status == 1
                unzip(filename);
%                 disp(['download of ', filename, ' successfully finished']);
                sFlag = 1;
            else
                set = set - 1; % try to access the latest file decrementally
                if set<1
                    sFlag=1;
                    
                end
            end
        end
end

% turn off stupid warning
warning('off', 'map:geotiff:undefinedGTModelTypeGeoKey');


% pass the eNDVI matrix.
if set>0
[A] = geotiffread(strcat(regFile, setStr, datestr(inDate,'yy'),'m.tif'));
R = worldfileread(strcat(regFile, setStr, datestr(inDate,'yy'),'m.tfw'),'geographic',size(A));
else
A=[];
R=[];
% disp(['No image before ',datestr(inDate),' in ',datestr(inDate,'yyyy')]);
end

warning('on', 'map:geotiff:undefinedGTModelTypeGeoKey');







