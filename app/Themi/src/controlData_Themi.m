
function controlData_Themi
% function to send out overview email.
%
% File:          controlData_Themi.m
%
% Created:        05/02/2015
%
% Last modified:  23/04/2015
%
% Author:         Sebastian Stoll, Jules enze, Tobias  (hydrosolutions ltd.)
%
% Purpose:        Script that assimilates the states and run the model.
%
% Description:   Script that assimilates the states and run the model.
%
%
% Copyright (C) 2015 hydrosolutions ltd. Zurich, Switzerland
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details. 


%clear all

tobiasTest = 0;

%Today's date
today=datenum(date);
%load setup file
load('setup.mat')

% if tobiasTest
%     setup.Path = '~/Dropbox (hydrosolutions)/iMoMoTOOLS/rrm/app';
%     setup.mPath = '~/Dropbox (hydrosolutions)/iMoMoTOOLS/rrm/app/Themi';
% end

cd(setup.Path)
cd('../src/dbase/')

%% To the JAVA thing
login = 'hydrosolution';
password = 'd8WmFgEGpZ';
url = 'https://157.26.64.27/HydrosolutionWS_V2/HydrosolutionWS';
namespace= 'http://webservice/';

javaclasspath('matlabSoapSecurity.jar');
secureSOAP = SecureSOAP(login, password);

%% Get daily values
yesterday=datestr(datenum(date)-1,'yyyy-mm-dd');
values = secureSOAP.callHTTPSSoap(url, namespace, 'getDateValues',{yesterday},{'date'}, {'{http://www.w3.org/2001/XMLSchema}string'});

% sites = [10,13,14,15,16,19,20,21,22,23,25,26,32,37,43,44,45,48,50,53,56,57,58,59,60,61,62,63,64,65,66,69];
% 
% for l=1:length(sites)
%     SID=find([values.siteID]==sites(l));
%     countS(l)=length(values(SID));
% end

% Filter values

id = find([values.siteID]==25);
count25 = length(values(id));

id = find([values.siteID]==58); 
count58 = length(values(id));

id = find([values.siteID]==55); 
count55 = length(values(id));

id = find([values.siteID]==59);
count59 = length(values(id));

id = find([values.siteID]==32); 
count32 = length(values(id));

id = find([values.siteID]==60); 
count60 = length(values(id));

id = find([values.siteID]==53); 
count53 = length(values(id));

id = find([values.siteID]==56); 
count56 = length(values(id));

id = find([values.siteID]==57); 
count57 = length(values(id));

id = find([values.siteID]==61); 
count61 = length(values(id));

id = find([values.siteID]==62); %Moivo Mollel - Selemani 0756541574
count62 = length(values(id));

id = find([values.siteID]==63);
count63 = length(values(id));

id = find([values.siteID]==64); % Mosses, Mwandeti furrow
count64 = length(values(id));

id = find([values.siteID]==65);
count65 = length(values(id));

id = find([values.siteID]==66); % Clemens, Urangini Furrow
count66 = length(values(id));

id = find([values.siteID]==22); % Clemens, Urangini Furrow
count22 = length(values(id));

id = find([values.siteID]==23); % Clemens, Urangini Furrow
count23 = length(values(id));

id = find([values.siteID]==37); % Clemens, Urangini Furrow
count37 = length(values(id));

id = find([values.siteID]==21); % Clemens, Urangini Furrow
count21 = length(values(id));

id = find([values.siteID]==50); % Clemens, Urangini Furrow
count50 = length(values(id));

id = find([values.siteID]==16); % Clemens, Urangini Furrow
count16 = length(values(id));

id = find([values.siteID]==43); % Clemens, Urangini Furrow
count43 = length(values(id));

id = find([values.siteID]==15); % Clemens, Urangini Furrow
count15 = length(values(id));

id = find([values.siteID]==44); % Clemens, Urangini Furrow
count44 = length(values(id));

id = find([values.siteID]==26); % Clemens, Urangini Furrow
count26 = length(values(id));

id = find([values.siteID]==45); % Clemens, Urangini Furrow
count45 = length(values(id));

id = find([values.siteID]==48); % Clemens, Urangini Furrow
count48 = length(values(id));

id = find([values.siteID]==14); % Clemens, Urangini Furrow
count14 = length(values(id));

id = find([values.siteID]==13); % Clemens, Urangini Furrow
count13 = length(values(id));

id = find([values.siteID]==10); % Bachuta, Moshi Weatherstation
count10 = length(values(id));

id = find([values.siteID]==35); % Bachuta, Moshi Soil Moisture
count35 = length(values(id));

id = find([values.siteID]==69); % Clemens, Urangini Furrow
count69 = length(values(id));


%% Load old data
cd(strcat(setup.mPath,'/data/processed/db/'))
delimiter = ',';
startRow = 2;
filename='report.csv';
formatSpec = '%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

for l=1:length(dataArray)
    data(:,l)=dataArray{l};
end

%Data
data{1,end+1}=datestr(today);
data{2,end}=count25;
data{3,end}=count58;
data{4,end}=count55;
data{5,end}=count59;
data{6,end}=count32;
data{7,end}=count60;
data{8,end}=count53;
data{9,end}=count56;
data{10,end}=count57;
data{11,end}=count61;
data{12,end}=count62;
data{13,end}=count63;
data{14,end}=count64;
data{15,end}=count65;
data{16,end}=count66;
data{17,end}=count22;
data{18,end}=count23;
data{19,end}=count37;
data{21,end}=count21;
data{22,end}=count50;
data{23,end}=count16;
data{24,end}=count43;
data{25,end}=count15;
data{26,end}=count44;
data{27,end}=count26;
data{28,end}=count45;
data{29,end}=count48;
data{30,end}=count14;
data{31,end}=count13;
data{32,end}=count10;
data{33,end}=count35;
data{34,end}=count69;

%% Write csv
% Details
data{1,1}='Site ID';data{1,2}='Person';data{1,3}='Telephone';data{1,4}='Location';

data{2,1}='25';data{2,2}='Lembros';data{2,3}='255765868402';data{2,4}='Olungura Siwandeti';

data{3,1}='58';data{3,2}='n/s';data{3,3}='n/s';data{3,4}='Nsanya Furrow';

data{4,1}='55';data{4,2}='n/s';data{4,3}='n/s';data{4,4}='Seliani River 1 at Nsanya';

data{5,1}='59';data{5,2}='n/s';data{5,3}='n/s';data{5,4}='Elakunoto Furrow';

data{6,1}='32';data{6,2}='Osambi';data{6,3}='255763469951';data{6,4}='Saitabau';

data{7,1}='60';data{7,2}='Osambi';data{7,3}='255763469951';data{7,4}='Olmotonyi Furrow';

data{8,1}='53';data{8,2}='Osambi';data{8,3}='255763469951';data{8,4}='Seliani River 2 at Josho';

data{9,1}='56';data{9,2}='Osambi';data{9,3}='255763469951';data{9,4}='Seliani River 3 at Kimunyaki';

data{10,1}='57';data{10,2}='Osambi';data{10,3}='255763469951';data{10,4}='Seliani RIver 4 at Lomunyaki';

data{11,1}='61';data{11,2}='Osambi';data{11,3}='255763469951';data{11,4}='Lomunyaki Furrow';

data{12,1}='62';data{12,2}='Osambi';data{12,3}='255763469951';data{12,4}='Ndeoya Furrow';

data{13,1}='63';data{13,2}='Osambi';data{13,3}='255763469951';data{13,4}='Kimunyaki Furrow';

data{14,1}='64';data{14,2}='Osambi';data{14,3}='255763469951';data{14,4}='Mwandeti Furrow at Seliani';

data{15,1}='65';data{15,2}='Osambi';data{15,3}='255763469951';data{15,4}='Olevolosi Furrow';

data{16,1}='66';data{16,2}='Osambi';data{16,3}='255763469951';data{16,4}='Oloirieni Furrow';

data{17,1}='22';data{17,2}='Osambi';data{17,3}='255763469951';data{17,4}='Olangit Olmringa/ Josho';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data{18,1}='23';data{18,2}='Osambi';data{18,3}='255763469951';data{18,4}='Saitabau/saumu';

data{19,1}='37';data{19,2}='Selemani';data{19,3}='255756541574';data{19,4}='Moivo Mollel';

data{21,1}='21';data{21,2}='Stephano';data{21,3}='255752291673';data{21,4}='Mfereji wa Wazee (olkajuloomoruack)';

data{22,1}='50';data{22,2}='Stephano';data{22,3}='255752291673';data{22,4}='Asa (shamba la mbegu)';

data{23,1}='16';data{23,2}='Thomas';data{23,3}='255759480365';data{23,4}='Kichangani Furrow';

data{24,1}='43';data{24,2}='Thomas';data{24,3}='255759480365';data{24,4}='Kigongoni Furrow';

data{25,1}='15';data{25,2}='Julius';data{25,3}='255768582513';data{25,4}='Mungushi Furrow';

data{26,1}='44';data{26,2}='Mosses';data{26,3}='255654106858';data{26,4}='Mwandeti Furrow';

data{27,1}='26';data{27,2}='Mussa';data{27,3}='255765826995';data{27,4}='Themi ya Simba Furrow';

data{28,1}='45';data{28,2}='Mussa';data{28,3}='255765826995';data{28,4}='Fili estate Furrow';

data{29,1}='48';data{29,2}='Clemens';data{29,3}='255752095454';data{29,4}='Urangini River';

data{30,1}='14';data{30,2}='Salehe';data{30,3}='255762162191';data{30,4}='Olokii School';

data{31,1}='13';data{31,2}='Bura Stanslavy';data{31,3}='255784946878';data{31,4}='Enyoito';

data{32,1}='10';data{32,2}='Bachuta';data{32,3}='255687112233';data{32,4}='Moshi-Weatherstation';

data{33,1}='35';data{33,2}='Bachuta';data{33,3}='255687112233';data{33,4}='PBWB Moshi - Soil Moisture';

data{34,1}='69';data{34,2}='David Charles';data{34,3}='255759916835';data{34,4}='TPC plantaion';


cd(strcat(setup.mPath,'/data/processed/db/'))
d=cell2table(data);
writetable(d,'report.csv');

%Settings
sender = 'imomowb@yahoo.com';
psswd = 'iMoMo567';

if ~tobiasTest
    recipient = {'rnaudascher@yahoo.de','siegfried@hydrosolutions.ch',...
        'tjkitomari@yahoo.com','bambabakary@yahoo.com',...
        'hoseasanga@gmail.com','alfayo.miseyeki@yahoo.com', 'patrice.mueller@he-arc.ch',...
        'felix.gaya-farre@hec.edu','marti@hydrosolutions.ch'};
else
    recipient = {'siegfried@hydrosolutions.ch', 'rnaudascher@yahoo.de'}; % testing
end

% YAHOO PREFERENCES
setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','smtp.mail.yahoo.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
% props.put('mail.smtp.starttls.enable','true');
%% 
% Message
count_site_ID=size(data,1)-1;
count_missing=0;
message=['The following sites did not send any data:' 10 ]
for i=2:size(data,1)
    switch data{i,1}
        case '25' %first site ID of class "Upstream" in the list above.
            message=[message 'Ngarenaro River:' 10];
        case '58'
            message=[message 'Seliani River:' 10];
        case '37'
            message=[message 'Themi River Upstream:' 10];
        case '21'
            message=[message 'Ngaramtoni River:' 10];
        case '16'
            message=[message 'Themi Downstream:' 10];
        case '14'
            message=[message 'Weather Stations / Mini-Cactus:' 10]; 
        case '69'
            message=[message 'TPC Plantation:' 10]; 
    end
    if data{i,end}==0
        message = [message '   ' 'Site ID' ' ' data{i,1} ': ' data{i,2} ' ' 'at' ' ' data{i,4} ' ' 'did not sent' 10];
        count_missing=count_missing+1;
    end
end
count_send=count_site_ID-count_missing;
message = ['Yesterday, ' values(1).dateTimeUTC(1:10) ', ' num2str(count_send) ' out of ' num2str(count_site_ID) ' sites submitted data' 10 message];

   
subject = 'Daily report on data sending activity';
attachment = {strcat('report.csv')};
sendmail(recipient, subject, message, attachment);  

end







