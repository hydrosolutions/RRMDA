

pfad='C:\Users\Sebastian\Dropbox (hydrosolutions)\iMoMoPHASE1\model\iMoMoTB_TS_NOV2014\app\themi'; %in future read this information from file
cd(pfad)


%% Do the JAVA thing
cd('../../src/dbase')

login = 'hydrosolution';
password = 'd8WmFgEGpZ';
url = 'https://157.26.64.27/HydrosolutionWS_V2/HydrosolutionWS';
namespace= 'http://webservice/';

javaclasspath('matlabSoapSecurity.jar');
secureSOAP = SecureSOAP(login, password);

%% Get daily values
values = secureSOAP.callHTTPSSoap(url, namespace, 'getDailyValues',{},{}, {});

%% Organize according to stations

out = organizeiMoMo(values,2);






% 
% 
% %Variable IDs
% vP=18; % Precipitation Eigentlich 59;
% vT=38; % Temperature
% vW=16; % Windspeed Eigentlich 54
% 
% 
% % Filter values
% idP=find([values.variableID]==vP);%Precipitation
% idT=find([values.variableID]==vT);%Temperature
% idW=find([values.variableID]==vW);%Windspeed
% id=[idP idT idW];
% data=values(id);
% 
% % Get SiteID and generate datenum
% for l=1:length(data)
%     sID(l)=data(l).siteID;
%     data(l).date=floor(datenum(data(l).dateTimeUTC,'yyyy-mm-dd HH:MM:SS'));
% end
% ID=unique(sID); % unique SiteID
% data = rmfield(data,'dateTimeUTC'); %Delete string date
% data=cell2mat(struct2cell(data)'); %Convert structure to matrix
% 
% % Get coordinates
% for l=1:length(ID)
%     LONLAT(l,1)=data(l,4);
%     LONLAT(l,2)=data(l,3);
% end
% 
% % Allocate values to each station and variable
% for l=1:length(ID)
%     id=find(data(:,5)==ID(l));
%     dummy=data(id,:);
%     id1=find(dummy(:,2)==vP);
%     dataS{l}.P=dummy(id1,:);
%     dataS{l}.dailyP=cell2mat(pivottable(num2cell(dataS{l}.P),[6],[],[1], @mean));
%     id2=find(dummy(:,2)==vT);
%     dataS{l}.T=dummy(id2,:);
%     dataS{l}.dailyT=cell2mat(pivottable(num2cell(dataS{l}.T),[6],[],[1], @mean));
%     dataS{l}.dailyMax=cell2mat(pivottable(num2cell(dataS{l}.T),[6],[],[1], @max));
%     dataS{l}.dailyMin=cell2mat(pivottable(num2cell(dataS{l}.T),[6],[],[1], @min));
%     id3=find(dummy(:,2)==vW);
%     dataS{l}.W=dummy(id3,:);
%     dataS{l}.dailyW=cell2mat(pivottable(num2cell(dataS{l}.W),[6],[],[1], @mean));
%        
% end
% 
% % Generate VARIABLES
% dateData=unique(data(:,6));
% YEARMODA=ones(length(dateData),length(ID)).*repmat(dateData,1,length(ID));
% dateN=zeros(length(dateData),length(ID));
% 
% for l=1:length(ID)
%     for k=1:length(dateData)
%         
%     %PRCP
%         if isempty(dataS{l}.dailyP)==1
%             PRCP(k,l)=NaN;
%             dateN(k,l)=dateN(k,l)+0;
%         else
%             flag=0;
%             for m=1:size(dataS{l}.dailyP,1)
%                 if dataS{l}.dailyP(m,1)==dateData(k)
%                     PRCP(k,l)=dataS{l}.dailyP(m,2);
%                     dateN(k,l)=dateN(k,l)+1;
%                     flag=1;
%                 else
%                     PRCP(k,l)=NaN;
%                     dateN(k,l)=dateN(k,l)+0;
%                 end
%                 if flag == 1
%                     break;
%                 end
%             end
%         end
%         
%      %TEMP
%         if isempty(dataS{l}.dailyT)==1
%             TEMP(k,l)=NaN;
%             dateN(k,l)=dateN(k,l)+0;
%         else
%             flag=0;
%             for m=1:size(dataS{l}.dailyT,1)
%                 if dataS{l}.dailyT(m,1)==dateData(k)
%                     TEMP(k,l)=dataS{l}.dailyT(m,2);
%                     dateN(k,l)=dateN(k,l)+1;
%                     flag=1;
%                 else
%                     TEMP(k,l)=NaN;
%                     dateN(k,l)=dateN(k,l)+0;
%                 end
%                 if flag == 1
%                     break;
%                 end
%             end
%         end
%         
%       %MAX
%         if isempty(dataS{l}.dailyMax)==1
%             MAX(k,l)=NaN;
%             dateN(k,l)=dateN(k,l)+0;
%         else
%             flag=0;
%             for m=1:size(dataS{l}.dailyMax,1)
%                 if dataS{l}.dailyMax(m,1)==dateData(k)
%                     MAX(k,l)=dataS{l}.dailyMax(m,2);
%                     dateN(k,l)=dateN(k,l)+1;
%                     flag=1;
%                 else
%                     MAX(k,l)=NaN;
%                     dateN(k,l)=dateN(k,l)+0;
%                 end
%                 if flag == 1
%                     break;
%                 end
%             end
%         end
%         
%        %MIN
%         if isempty(dataS{l}.dailyMin)==1
%             MIN(k,l)=NaN;
%             dateN(k,l)=dateN(k,l)+0;
%         else
%             flag=0;
%             for m=1:size(dataS{l}.dailyMin,1)
%                 if dataS{l}.dailyMin(m,1)==dateData(k)
%                     MIN(k,l)=dataS{l}.dailyMin(m,2);
%                     dateN(k,l)=dateN(k,l)+1;
%                     flag=1;
%                 else
%                     MIN(k,l)=NaN;
%                     dateN(k,l)=dateN(k,l)+0;
%                 end
%                 if flag == 1
%                     break;
%                 end
%             end
%         end
%         
%         %WDSP
%         if isempty(dataS{l}.dailyW)==1
%             WDSP(k,l)=NaN;
%             dateN(k,l)=dateN(k,l)+0;
%         else
%             flag=0;
%             for m=1:size(dataS{l}.dailyW,1)
%                 if dataS{l}.dailyW(m,1)==dateData(k)
%                     WDSP(k,l)=dataS{l}.dailyW(m,2);
%                     dateN(k,l)=dateN(k,l)+1;
%                     flag=1;
%                 else
%                     WDSP(k,l)=NaN;
%                     dateN(k,l)=dateN(k,l)+0;
%                 end
%                 if flag == 1
%                     break;
%                 end
%             end
%         end  
%         
%     end
% end
% YEARMODA(dateN<=0)=NaN;
% 










































