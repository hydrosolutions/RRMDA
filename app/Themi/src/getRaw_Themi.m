function getRaw_Themi
% Script to download different raw data from the internet.
%
% File:           getRaw_Themi.m
%
% Created:        13/12/2012
%
% Last modified:  23/04/2015 TS: Stabilization and error catching
%
% Author:         Sebastian Stoll (hydrosolutions ltd.)
%                 Tobias Siegfried
%
% Purpose:        Download different raw data and save them.
%
% Description:    Download different raw data and save them.
%
%
% Copyright (C) 2015 hydrosolutions ltd. Zurich, Switzerland
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.

try
    
    %% Spatial and temporal definitions
    tic
    clear all
    warning off
    
    %Today's date
    today = datenum(date);
    
    % Open log file for writing.
    try
        fid = fopen('getRaw_Themi.log','a+');
    catch ME
        disp('Error opening log file for writing.')
        rethrow(ME)
    end
    c = clock;
    doublefprintf(fid,'-------------------------------------------------\n');
    doublefprintf(fid,' %02d.%02d.%04d %02d:%02d : Starting getRaw_Themi\n',c(3),c(2),c(1),c(4),c(5));
    doublefprintf(fid,'-------------------------------------------------\n');
    
    doublefprintf(fid,'Setting up . . . ')
    %load setup file
    load('setup.mat')
    
    %bounding box
    bblon = [setup.bBox(1) setup.bBox(3)];
    bblat = [setup.bBox(2) setup.bBox(4)];
    ndays = 6; %number of days for download of GDAS/GFS (Maximum = 6)-1
    pfad = setup.mPath;
    region = setup.ndvi; %in future read this information from file
    cd(pfad)
    doublefprintf(fid,'. . . done.\n');
    
    %% 1. FEWS PRECIPITATION in mm/d
    doublefprintf(fid,'Getting precipitation data . . . ');
    cd('./data/raw')
    
    % read new data
    [new.P,new.R_P,new.timeP] = eval(strcat('get',setup.P,'(today,bblat,bblon)'));
    if isempty(new.P)
        doublefprintf(fid,'ERROR: Failed in call get%s(%d,[%d,%d],[%d,%d]).\nContinguing.',setup.P,today,bblat,bblon);
    end
    
    % check if there is already a file
    if exist(fullfile(cd, strcat(setup.P,'.mat')), 'file')==0
        
        %if there is no old FEWS file, save the new data in the following variables
        P = new.P;
        R_P = new.R_P;
        timeP = new.timeP;
        
    else
        %if there is an old FEWS file, open it
        load(strcat(setup.P,'.mat'));
        
        %Update the data
        if max(timeP,[],2) ~= max(new.timeP,[],2) % if there is no new data, save the old data
            nNew = max(new.timeP,[],2)-max(timeP,[],2); % number of days of new data
            if nNew > length(new.timeP) % Insert NaNs
              nMiss = nNew - length(new.timeP);
              % fill in missing values
              new.timeP = (new.timeP(1)-nMiss:new.timeP(end))';
              P = cat(3,P,NaN(size(P,1),size(P,2),nMiss));
              nNew = nNew - nMiss;
              timeP(end+1:end+nNew) = new.timeP(end-nNew+1:end); % add new dates
              P(:,:,end+1:end+nNew) = new.P(:,:,end-nNew+1:end); % add new temp data
            else
              timeP(end+1:end+nNew) = new.timeP(end-nNew+1:end); % add new dates
              P(:,:,end+1:end+nNew) = new.P(:,:,end-nNew+1:end); % add new precipitation data
            end
        end
    end
    
    %Save updated data
    save(strcat(setup.P,'.mat'),'P','R_P','timeP');
    %disp('= FEWS P loaded')
    
    %Clear variables not needed anymore
    clearvars P R_P timeP new
    doublefprintf(fid,'. . . done.\n');
    
    %% 2. GDAS Temperatures
    doublefprintf(fid,'Getting temperature data . . . ');
    
    % read new data
    [new.T,new.maxT,new.minT,new.R_T,new.timeT] = eval(strcat('get',setup.T,'(today-ndays,today,bblat,bblon)'));
    if isempty(new.T)
        doublefprintf(fid,'ERROR: Failed in call to get%s(%d,%d,[%d,%d],[%d,%d]).\nContinuing.\n',setup.T,today-ndays,today,bblat,bblon);
    end
    
    % check if there is already a file
    if exist(fullfile(cd, strcat(setup.T,'.mat')), 'file')==0
        
        %if there is no old GDAS file, save the new data in the following variables
        T = new.T;
        maxT = new.maxT;
        minT = new.minT;
        R_T = new.R_T;
        timeT = new.timeT;
        
    else
        %if there is an old GDAS file, open it
        load(strcat(setup.T,'.mat'));
        
        %Update the data
        if max(timeT,[],2) ~= max(new.timeT,[],1) % if there is no new data, save the old data
            nNew = max(new.timeT,[],1)-max(timeT,[],2); % number of new data
            % check if data gap exists and correct! (TS: 22/04/2015)
            if nNew <= ndays
                
                timeT(end+1:end+nNew) = new.timeT(end-nNew+1:end); % add new dates
                T(:,:,end+1:end+nNew) = new.T(:,:,end-nNew+1:end); % add new temp data
                maxT(:,:,end+1:end+nNew) = new.maxT(:,:,end-nNew+1:end); % add new temp data
                minT(:,:,end+1:end+nNew) = new.minT(:,:,end-nNew+1:end); % add new temp data
                
            else % more missing days than GDAS data - do NaN tagging
                
                nMiss = nNew - ndays;
                % fill in missing values
                new.timeT = (new.timeT(1)-nMiss:new.timeT(end))';
                T = cat(3,T,NaN(size(T,1),size(T,2),nMiss));
                maxT = cat(3,maxT,NaN(size(T,1),size(T,2),nMiss));
                minT = cat(3,minT,NaN(size(T,1),size(T,2),nMiss));
                nNew = ndays;
                timeT(end+1:end+nNew) = new.timeT(end-nNew+1:end); % add new dates
                T(:,:,end+1:end+nNew) = new.T(:,:,end-nNew+1:end); % add new temp data
                maxT(:,:,end+1:end+nNew) = new.maxT(:,:,end-nNew+1:end); % add new max temp data
                minT(:,:,end+1:end+nNew) = new.minT(:,:,end-nNew+1:end); % add new min temp data
                
            end
        end
    end
    
    %Save updated data
    save(strcat(setup.T,'.mat'),'T','maxT','minT','R_T','timeT');
    %disp('= GDAS T loaded')
    
    %Clear variables not needed anymore
    clearvars T maxT minT timeT new nNew R_T
    
    doublefprintf(fid,'. . . done.\n');
    
    %% 3. GFS Precipitation and temperature forecasts
    doublefprintf(fid,'Getting forecast data . . . ');
    %minimum lat/lon difference =1 degree
    GFSbblat(1) = floor(bblat(1));
    GFSbblat(2) = ceil(bblat(2));
    GFSbblon(1) = floor(bblon(1));
    GFSbblon(2) = ceil(bblon(2));
    
    % read new data
    [new] = eval(strcat('get',setup.F,'(today-ndays,today,GFSbblat,GFSbblon)'));
    if isempty(new)
        doublefprintf(fid,'ERROR: Failed in call to get%s(%d,%d,[%d,%d],[%d,%d]).\nContinuing.\n',setup.F,today-ndays,today,GFSbblat,GFSbblon);
    end
    
    % check if there is already a file
    if exist(fullfile(cd, strcat(setup.F,'.mat')), 'file')==0
        
        %if there is no old GFS file, save the new data in the following variables
        F = new;
        
    else
        
        %if there is an old GFS file, open it
        load(strcat(setup.F,'.mat'));
        
        if F{end}.timeF(1) ~= new{end}.timeF(1);
            
            nNew = new{end}.timeF(1)-F{end}.timeF(1); % number of new data
            
            if nNew <= ndays
                
                % Simply update the data
                F(end+1:end+nNew) = new(end-nNew+1:end);
                
            else
                
                % NaN-tag and update the data
                nMiss = nNew - ndays;
                
                newMissing = new(1);
                
                tempTimeF = newMissing{1}.timeF;
                
                for idx = 1 : nMiss
                    newMissing{idx}.FT = NaN * newMissing{1}.FT;
                    newMissing{idx}.FmaxT = NaN * newMissing{1}.FmaxT;
                    newMissing{idx}.FminT = NaN * newMissing{1}.FminT;
                    newMissing{idx}.FP = NaN * newMissing{1}.FP;
                    newMissing{idx}.timeF = tempTimeF - ((nMiss - idx) + 1);
                end
                
                new = [newMissing new];
                
                F(end+1:end+nNew) = new(end-nNew+1:end);
                
            end
            
        end
    end
    
    %Georeference of forecast
    R_F = F{end}.R_F;
    
    %Save updated data
    save(strcat(setup.F,'.mat'),'F');
    %disp('= GFS P/T loaded')
    
    %Clear variables not needed anymore
    clearvars F new nNew
    doublefprintf(fid,'. . . done.\n');
    
    %% 4. Get WMO station data
    doublefprintf(fid,'Getting WMO station data . . . ');
    % create temporary directory and change dir
    tempDir = 'temp'; % temporary directory for data download
    removeTemp = 1;   % remove temporary directory: yes = 1, no = 0;
    mkdir(tempDir);
    cd(tempDir);
    
    % get weather station list (global)
    cd('../../../../../src/data')
    load st.mat;
    load glocloud.mat;
    
    cd(pfad)
    cd('./data/raw/temp')
    
    % extent for weather stations to be downloaded
    ll = [GFSbblon(1) GFSbblat(1)];           % lower left corner
    ur= [GFSbblon(2) GFSbblat(2)];             % upper right corner
    
    % read data
    new = getStations(ll, ur, st, today, 20);
    
    %flag NaNs in PRCP
    P=new(1).PRCP;
    P(P==99.9900)=NaN;
    new(1).PRCP=P;
    
    % count and delete last days without date
    nsum = sum(isnan(new.YEARMODA),2);
    isum = find(nsum==size(new.YEARMODA,2));
    new.YEARMODA(isum,:) = [];
    new.TEMP(isum,:) = [];
    new.DEWP(isum,:) = [];
    new.WDSP(isum,:) = [];
    new.MAX(isum,:) = [];
    new.MIN(isum,:) = [];
    new.PRCP(isum,:) = [];
    new.FRSHTT(isum,:) = [];
    
    cd('..');
    if removeTemp == 1
        rmdir(tempDir, 's')
    end
    
    % calculate potET and set units right
    new.nN = cCloud(new.LONLAT(:,1),new.LONLAT(:,2), glocloud); % LON, LAT, cloud
    new.pET = etStations(new);
    new.TEMP = (new.TEMP - 32).*5./9; %Fahrenheit to Celcius
    new.MAX = (new.MAX - 32).*5./9; %Fahrenheit to Celcius
    new.MIN = (new.MIN - 32).*5./9; %Fahrenheit to Celcius
    new.PRCP = new.PRCP.*25.4; %Inch to Millimeters
    
    % check if there is already a file
    if exist(fullfile(cd, 'WMO.mat'), 'file')==0
        
        %if there is no old GFS file, save the new data in the following variables
        WMO = new;
        
    else
        %if there is an old GFS file, open it
        load('WMO.mat');
        
        %find row with least NaN
        a = WMO.YEARMODA;
        a(isnan(a)==0) = 0;
        a(isnan(a)==1) = 1;
        a = sum(a,1);
        [v ID] = min(a);
        
        % if there is new data update the data
        if isnan(new.YEARMODA(end,ID))==0
            if WMO.YEARMODA(end,ID) ~= new.YEARMODA(end,ID) ; % if there is no new data, save the old data
                % number of new data
                nNew=datenum(num2str(new.YEARMODA(end,ID)),'yyyymmdd')-datenum(num2str(WMO.YEARMODA(end,ID)),'yyyymmdd');
                if nNew > size(new.YEARMODA,1)  % Patch with NaNs
                  nMiss = nNew - size(new.YEARMODA,1);
                  for l = 1:nMiss
                    WMO.YEARMODA(end+l,:) = NaN .* new.YEARMODA(end,:);
                    WMO.TEMP(end+l,:) = NaN .* new.TEMP(end,:);
                    WMO.DEWP(end+l,:) = NaN .* new.DEWP(end,:);
                    WMO.WDSP(end+l,:) = new.WDSP(end,:);
                    WMO.MAX(end+l,:) = new.MAX(end,:);
                    WMO.MIN(end+l,:) = new.MIN(end,:);
                    WMO.PRCP(end+l,:) = new.PRCP(end,:);
                    WMO.FRSHTT(end+l,:) = new.FRSHTT(end,:);
                  end
                  nNew = nNew - nMiss;
                  WMO.YEARMODA(end+1:end+nNew,:) = new.YEARMODA(end-nNew+1:end,:);
                  WMO.TEMP(end+1:end+nNew,:) = new.TEMP(end-nNew+1:end,:);
                  WMO.DEWP(end+1:end+nNew,:) = new.DEWP(end-nNew+1:end,:);
                  WMO.WDSP(end+1:end+nNew,:) = new.WDSP(end-nNew+1:end,:);
                  WMO.MAX(end+1:end+nNew,:) = new.MAX(end-nNew+1:end,:);
                  WMO.MIN(end+1:end+nNew,:) = new.MIN(end-nNew+1:end,:);
                  WMO.PRCP(end+1:end+nNew,:) = new.PRCP(end-nNew+1:end,:);
                  WMO.FRSHTT(end+1:end+nNew,:) = new.FRSHTT(end-nNew+1:end,:);
                else
                  %add new data
                  WMO.YEARMODA(end+1:end+nNew,:) = new.YEARMODA(end-nNew+1:end,:);
                  WMO.TEMP(end+1:end+nNew,:) = new.TEMP(end-nNew+1:end,:);
                  WMO.DEWP(end+1:end+nNew,:) = new.DEWP(end-nNew+1:end,:);
                  WMO.WDSP(end+1:end+nNew,:) = new.WDSP(end-nNew+1:end,:);
                  WMO.MAX(end+1:end+nNew,:) = new.MAX(end-nNew+1:end,:);
                  WMO.MIN(end+1:end+nNew,:) = new.MIN(end-nNew+1:end,:);
                  WMO.PRCP(end+1:end+nNew,:) = new.PRCP(end-nNew+1:end,:);
                  WMO.FRSHTT(end+1:end+nNew,:) = new.FRSHTT(end-nNew+1:end,:);
                end
            end
        end
    end
    
    
    %Save updated data
    save('WMO.mat','WMO');
    %disp('= WMO Stations loaded')
    
    %Clear variables not needed anymore
    clearvars P WMO new nNew a
    doublefprintf(fid,'. . . done.\n');
    
    
    %% 5. GET eNDVI
    % reads the newest eNDVI image and saves it with the date as filename
    doublefprintf(fid,'Getting eNDVI data . . . ');
    
    cd('./ndvi')
    
    if strcmp(setup.ndvi,'N/A')==0
        
        % get NDVI date reference
        load ref.mat;
        ref(:,1)=str2double(datestr(today,'yyyy')); %set current year
        ref=datenum(ref);% serial date
        
        % create temporary directory and change dir
        tempDir = 'temp'; % temporary directory for data download
        removeTemp = 1;   % remove temporary directory: yes = 1, no = 0;
        mkdir(tempDir);
        cd(tempDir);
        
        %reads the latest eNDVI image and cut it to GDAS Dimensions
        try
            [NDVI,set,R_NDVI] = getNDVI(today,region);
            
            try
                [NDVI,R_NDVI] = maptrims(NDVI,R_NDVI,R_F.LatitudeLimits,R_F.LongitudeLimits);
            catch matlabError2
                doublefprintf(fid,'\nERROR: Failed in call to maptrims.\n')
                errorMessage = getReport(matlabError2);
                doublefprintf(fid,'%s\n',errorMessage);
                %rethrow(ME)
            end
            
            % remove temporary folder
            cd('..');
            if removeTemp == 1
                rmdir(tempDir, 's')
            end
        
            if isempty(NDVI)==0
                % date of image
                timeNDVI=ref(set);
            
                % Save new NDVI image
                save(strcat('NDVI_',num2str(timeNDVI),'.mat'),'NDVI','timeNDVI','R_NDVI');
                %disp('= eNDVI loaded')
            end
        
            clearvars NDVI set removeTemp tempDir ref timeNDVI R_NDVI R_T
        
        catch matlabError
            cd('..');
            doublefprintf(fid,'\nERROR: Failed in call to getNDVI(%d,%s).\n',today,region);
            errorMessage = getReport(matlabError);
            doublefprintf(fid,'%s\n',errorMessage);
            %rethrow(ME)
            
            
        end
        
        
    end
    doublefprintf(fid,'. . . done.\n');
    
    %% 6. Get iMoMo Station data
    doublefprintf(fid,'Getting iMoMo station data. . . ');
    
    % Do the JAVA thing
    cd('../../../../../src/dbase')
    
    login = 'hydrosolution';
    password = 'd8WmFgEGpZ';
    url = 'https://157.26.64.27/HydrosolutionWS_V2/HydrosolutionWS';
    namespace= 'http://webservice/';
    yesterday=datestr(datenum(date)-1,'yyyy-mm-dd');
    
    javaclasspath('matlabSoapSecurity.jar');
    secureSOAP = SecureSOAP(login, password);
    
    % Get yesterdays values
    values = secureSOAP.callHTTPSSoap(url, namespace, 'getDateValues',{yesterday},{'date'}, {'{http://www.w3.org/2001/XMLSchema}string'});
    
    %% Save csv files
    if ischar(values)==0
        %delete forecast entries
        id1=find([values.variableID]>9);
        id2=find([values.variableID]<28);
        id3=intersect(id1,id2);
        values(id3)=[];
    else
        values=[];
    end
    
    % Get site details
    for l=1:length(values)
        siteD{l} = secureSOAP.callHTTPSSoap(url, namespace, 'getSiteDetails',{values(l).siteID},{'siteID'}, {'{http://www.w3.org/2001/XMLSchema}int'});
        values(l).siteName=siteD{l}.siteName;
    end
    
    % Get user name
    load('user.mat')
    uID=cell2mat(user(:,1));
    uN=user(:,2);
    for l=1:length(values)
        name{l}=(uN(uID==values(l).userID));
        values(l).userName=name{l};
    end
    
    %save data and write csv
    cd(strcat(setup.mPath,'/data/raw/imomo/'))
    
    if isempty(values)==0
        DL(2:length(values)+1,:)=struct2cell(values)';
        %empty cells
        ec=cellfun(@isempty,DL);
        DL(ec)={'N/A'};
        DL(1,1)={'Data Value'};
        DL(1,2)={'Variable ID'};
        DL(1,3)={'Date'};
        DL(1,4)={'Latitude'};
        DL(1,5)={'Longitude'};
        DL(1,6)={'Site ID'};
        DL(1,7)={'User ID'};
        DL(1,8)={'Site Name'};
        DL(1,9)={'User Nickname'};
        f=cell2table(DL);
        writetable(f,strcat(num2str(today),'DLoad.csv'))
    end
    
    
    %% Find and save Discharge data
    data=values;
    cd(pfad)
    cd('data/raw/imomo')
    save(strcat(num2str(today),'DLoad'),'data');
    
    % Meteo Station list
    stations=[10 13 14];
    for l=1:length(stations)
        detailS(l) = secureSOAP.callHTTPSSoap(url, namespace, 'getSiteDetails',{stations(l)},{'siteID'}, {'{http://www.w3.org/2001/XMLSchema}int'});
    end
    
    % Organize Meteo data according to stations
    stations=[10 13 14];
    for l=1:length(stations)
        out{l}=organizeiMoMo(data,stations(l));
    end
    
    cd(pfad)
    cd('data/raw')
    % check if there is already a file
    if exist(fullfile(cd, 'iMoMo.mat'), 'file')==0
        
        %if there is no old iMoMo file, save the new data in the following variables
        for l=1:length(stations)
            iMoMo.ID(l,1)=stations(l);
            iMoMo.LONLAT(l,:)=[detailS(l).longitude detailS(l).latitude];
            iMoMo.PRCP(:,l)=out{l}.PRCP;
            iMoMo.TEMP(:,l)=out{l}.TEMP;
            iMoMo.MAX(:,l)=out{l}.MAX;
            iMoMo.MIN(:,l)=out{l}.MIN;
            iMoMo.WDSP(:,l)=out{l}.WDSP;
            iMoMo.YEARMODA(:,l)=out{l}.YEARMODA;
        end
        
        save iMoMo.mat iMoMo
        
    else
        %if there is an old iMoMo file, open it
        load('iMoMo.mat');
        nextD=size(iMoMo.PRCP,1)+1;
        numSt=size(iMoMo.ID,1);
        
        %check if it is the same number of stations
        if numSt==length(stations)
            for l=1:length(stations)
                %Fill in new data
                iMoMo.PRCP(nextD,l)=out{l}.PRCP(end);
                iMoMo.TEMP(nextD,l)=out{l}.TEMP(end);
                iMoMo.MAX(nextD,l)=out{l}.MAX(end);
                iMoMo.MIN(nextD,l)=out{l}.MIN(end);
                iMoMo.WDSP(nextD,l)=out{l}.WDSP(end);
                iMoMo.YEARMODA(nextD,l)=out{l}.YEARMODA(end);
            end
            save iMoMo.mat iMoMo
        else
            for l=1:length(stations)
                %Fill in new data
                iMoMo.PRCP(nextD,l)=out{l}.PRCP;
                iMoMo.TEMP(nextD,l)=out{l}.TEMP;
                iMoMo.MAX(nextD,l)=out{l}.MAX;
                iMoMo.MIN(nextD,l)=out{l}.MIN;
                iMoMo.WDSP(nextD,l)=out{l}.WDSP;
                iMoMo.YEARMODA(nextD,l)=out{l}.YEARMODA;
                iMoMo.ID(l,1)=stations(l);
                iMoMo.LONLAT(l,:)=[detailS(l).longitude detailS(l).latitude];
            end
            %Fill in NaN for new stations
            nNew=length(stations)-numSt; %new number of stations
            
            iMoMo.PRCP(1:nextD-1,end-(nNew-1):end)=NaN;
            iMoMo.TEMP(1:nextD-1,end-(nNew-1):end)=NaN;
            iMoMo.MAX(1:nextD-1,end-(nNew-1):end)=NaN;
            iMoMo.MIN(1:nextD-1,end-(nNew-1):end)=NaN;
            iMoMo.WDSP(1:nextD-1,end-(nNew-1):end)=NaN;
            iMoMo.YEARMODA(1:nextD-1,end-(nNew-1):end)=NaN;
            save iMoMo.mat iMoMo
        end
    end
    
    %disp('= iMoMo Data loaded')
    doublefprintf(fid,'. . . done.\n');
    
    
    %% Close all and display duration
    c = clock;
    doublefprintf(fid,'-------------------------------------------------\n');
    doublefprintf(fid,' %02d.%02d.%04d %02d:%02d : getRaw_Themi successfully finished!\n',c(3),c(2),c(1),c(4),c(5));
    doublefprintf(fid,'-------------------------------------------------\n');
    
    close all
    disp('---')
    disp([datestr(today) ': getRaw_Themi successfully finished!'])
    disp('---')
    
    toc
    
    
catch ME
    
    disp('---')
    disp('time')
    nowT = datevec(now);
    nowT = nowT(1:5)
    disp('---')
    
    fid = fopen('getRaw_Themi.log','a+');
    report = getReport(ME);
    doublefprintf(fid,'# %2d.%2d.%4d %2d:%2d\n',nowT(3),nowT(2),nowT(1),nowT(4),nowT(5));
    doublefprintf(fid,'%s',report);
    fclose(fid);
    
    disp('Error occured in getRaw_Themi.m! ERROR SPECIFICS:')
    rethrow(ME)
    disp('---')
    disp('Next try in 3 hours!')
    disp('---')
    
    
    
end

end % getRaw_Themi



%% Local functions

function doublefprintf(fid,message,varargin)
%% function doubpefprintf prints to file handle and to console.
%
% function doublefprintf(fid,message,varargin)
%
% @param fid int file handle to a writable file.
% @param message string with message to print.
% @param varargin arguments to be passed to the message.
%

if nargin < 3
    fprintf(fid,message);
    fprintf(message);
elseif nargin >= 3
    fprintf(fid,message,varargin{:});
    fprintf(message,varargin{:});
end

end























