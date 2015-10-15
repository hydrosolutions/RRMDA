function [data1] = getObservation(type,day,nC,paths)
%% function data1 = getObservation(type, day, nC, paths)
%
% Reads in the current iMoMo data and extracts discharge measurements for
% each subcatchment. 
%
% @arg type (string) indicating type of observed data. E.g. 'Q'. 
% @arg day  (int) current time step to compute.
% @arg nC   (int) number of sub-catchments in the model area. 
% @arg paths (struct) with paths to relevant directories.
% @return data1 (double array) with dimension [1,nC] containing the newest 
%           observation for each subcatchment.
%
% 

%Find file
cd(strcat(paths.data,'raw/imomo'))

if exist(strcat(num2str(day)), 'file')
   load(strcat(num2str(day),'DLoad'))

    %Find typeID
    switch type
        case 'Q'
            typeID=25;
    end


    % Get Values
    if length(data)>0 %is there is a new observation

        %Resrict to data type
        tID=find([data.variableID]==typeID);
        typeData=values(tID);

        for C=1:nC % for each subcatchment
            sID=find([typeData.siteID]==C);
            cData{C}=typeData(sID);

            if isempty(cData{C}) %if there is no observation for that subcatchment
                data1(C)=NaN;

            else % if there is a observation for that subcatchment get the latest observations.
                days=datenum({cData{C}.dateTimeUTC},'yyyy-mm-dd HH:MM:SS');
                [dum idD] = max(days);
                dum1=cData{C};
                data1(C)=dum1(idD).dataValue;
            end

        end

    else %if there is no observation
        data1(1:nC)=NaN;
    end
else
    data1(1:nC)=NaN;

end
        
    
    

  
