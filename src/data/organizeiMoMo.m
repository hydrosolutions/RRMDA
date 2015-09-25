


function out = organizeiMoMo(values,stationID)


%Variable IDs
vP=59; % Precipitation Eigentlich;
vT=38; % Temperature
vW=54; % Windspeed Eigentlich


% Filter values
idP=find([values.variableID]==vP);%Precipitation
idT=find([values.variableID]==vT);%Temperature
idW=find([values.variableID]==vW);%Windspeed
id=[idP idT idW];
data=values(id);


if isempty(data)==0
    

    %Select StationID
    idS=find([data.siteID]==stationID);%Windspeed
    data=data(idS);

    if isempty(data)==0

        % Generate datenum
        for l=1:size(data,1)
            data(l).date=floor(datenum(data(l).dateTimeUTC,'yyyy-mm-dd HH:MM:SS'));
        end
        
        data = rmfield(data,'dateTimeUTC'); %Delete string date
        data = rmfield(data,'userName'); %Delete string date
        data = rmfield(data,'siteName'); %Delete string date
        data=cell2mat(struct2cell(data)'); %Convert structure to matrix


        % Get coordinates
        LONLAT(1)=data(1,4);
        LONLAT(2)=data(1,3);

        % Allocate values to each station and variable
        id1=find(data(:,2)==vP);
        dataS.P=data(id1,:);
        dataS.dailyP=cell2mat(pivottable(num2cell(dataS.P),[6],[],[1], @mean));

        id2=find(data(:,2)==vT);
        dataS.T=data(id2,:);
        dataS.dailyT=cell2mat(pivottable(num2cell(dataS.T),[6],[],[1], @mean));
        dataS.dailyMax=cell2mat(pivottable(num2cell(dataS.T),[6],[],[1], @max));
        dataS.dailyMin=cell2mat(pivottable(num2cell(dataS.T),[6],[],[1], @min));

        id3=find(data(:,2)==vW);
        dataS.W=data(id3,:);
        dataS.dailyW=cell2mat(pivottable(num2cell(dataS.W),[6],[],[1], @mean));

        % Generate VARIABLES
        dateData=unique(data(:,6));
        YEARMODA=ones(length(dateData),1).*repmat(dateData,1,1);
        dateN=zeros(length(dateData),1);

        for k=1:length(dateData)

            %PRCP
                if isempty(dataS.dailyP)==1
                    PRCP(k)=NaN;
                    dateN(k)=dateN(k)+0;
                else
                    flag=0;
                    for m=1:size(dataS.dailyP,1)
                        if dataS.dailyP(m,1)==dateData(k)
                            PRCP(k)=dataS.dailyP(m,2);
                            dateN(k)=dateN(k)+1;
                            flag=1;
                        else
                            PRCP(k)=NaN;
                            dateN(k)=dateN(k)+0;
                        end
                        if flag == 1
                            break;
                        end
                    end
                end

             %TEMP
                if isempty(dataS.dailyT)==1
                    TEMP(k)=NaN;
                    dateN(k)=dateN(k)+0;
                else
                    flag=0;
                    for m=1:size(dataS.dailyT,1)
                        if dataS.dailyT(m,1)==dateData(k)
                            TEMP(k)=dataS.dailyT(m,2);
                            dateN(k)=dateN(k)+1;
                            flag=1;
                        else
                            TEMP(k)=NaN;
                            dateN(k)=dateN(k)+0;
                        end
                        if flag == 1
                            break;
                        end
                    end
                end

              %MAX
                if isempty(dataS.dailyMax)==1
                    MAX(k)=NaN;
                    dateN(k)=dateN(k)+0;
                else
                    flag=0;
                    for m=1:size(dataS.dailyMax,1)
                        if dataS.dailyMax(m,1)==dateData(k)
                            MAX(k)=dataS.dailyMax(m,2);
                            dateN(k)=dateN(k)+1;
                            flag=1;
                        else
                            MAX(k)=NaN;
                            dateN(k)=dateN(k)+0;
                        end
                        if flag == 1
                            break;
                        end
                    end
                end

               %MIN
                if isempty(dataS.dailyMin)==1
                    MIN(k)=NaN;
                    dateN(k)=dateN(k)+0;
                else
                    flag=0;
                    for m=1:size(dataS.dailyMin,1)
                        if dataS.dailyMin(m,1)==dateData(k)
                            MIN(k)=dataS.dailyMin(m,2);
                            dateN(k)=dateN(k)+1;
                            flag=1;
                        else
                            MIN(k)=NaN;
                            dateN(k)=dateN(k)+0;
                        end
                        if flag == 1
                            break;
                        end
                    end
                end

                %WDSP
                if isempty(dataS.dailyW)==1
                    WDSP(k)=NaN;
                    dateN(k)=dateN(k)+0;
                else
                    flag=0;
                    for m=1:size(dataS.dailyW,1)
                        if dataS.dailyW(m,1)==dateData(k)
                            WDSP(k)=dataS.dailyW(m,2);
                            dateN(k)=dateN(k)+1;
                            flag=1;
                        else
                            WDSP(k)=NaN;
                            dateN(k)=dateN(k)+0;
                        end
                        if flag == 1
                            break;
                        end
                    end
                end  


        end
        YEARMODA(dateN<=0)=NaN;


        out.YEARMODA=YEARMODA;
        out.MAX=MAX;
        out.MIN=MIN;
        out.PRCP=PRCP;
        out.TEMP=TEMP;
        out.WDSP=WDSP;
        out.ID=stationID;
        out.LONLAT=LONLAT;
    else
        out.YEARMODA=NaN;
        out.MAX=NaN;
        out.MIN=NaN;
        out.PRCP=NaN;
        out.TEMP=NaN;
        out.WDSP=NaN;
        out.ID=NaN;
        out.LONLAT=NaN;
    end

else
    out.YEARMODA=NaN;
    out.MAX=NaN;
    out.MIN=NaN;
    out.PRCP=NaN;
    out.TEMP=NaN;
    out.WDSP=NaN;
    out.ID=NaN;
    out.LONLAT=NaN;
end

   
end











































