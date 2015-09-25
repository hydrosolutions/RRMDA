function [rasterC, slope, Vi] = biasCorrection(raster, rasterR, timeT, statData, var2BC, dataSource)
% Bias correction of remotely-sensed data based on station data
%
% @raster           - raster data to be bias corrected
% @statData         - structure containing station data (iMoMo Met Stations and WMO stations
% @corrType         - Temperature (additive) or Precipitation (multiplicative)
% @return rasterC   - corrected raster data
%
% Usage:
%                   sVal = biasCorrection(raster, statData, corrType)
%
% File:             biasCorrection.m
%
% Created:          07/01/2015
%
% Last modified:    07/01/2015
%
% Author:           Tobias Siegfried (hydrosolutions ltd.)
%
% Purpose:          Bias correction of remotely-sensed data based on station
%                   data (based on inverse distance correction).
%
% Description:      Ditto
%
% Revisions:        NA, Note: there should be also an option to include
%                   station data from other sources (apart from iMoMo and
%                   WMO), e.g. from local met services
%
% Copyright:        2015, hydrosolutions ltd.
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.

%% BODY

debugging = 0;

if strcmp(dataSource, 'GFS')

    % Identify raster data set
    if strcmp(var2BC, 'TEMP')
        varF='FT';
    elseif strcmp(var2BC, 'MIN')
        varF='FminT';
    elseif strcmp(var2BC, 'MAX')
        varF='FmaxT';
    elseif strcmp(var2BC, 'PRCP')
        varF='FP';
    end
   
    sE=size(eval(strcat('raster{1}.',varF)),4);
    sD=size(eval(strcat('raster{1}.',varF)),3);
    rasterC=raster;
    
%for each ensemble
    for k=1:sE 
        % Generate time vector
        for l=1:timeT
                if isfield(raster{l},varF)==1
                    timeF(l)=raster{l}.timeF(1);
                else
                    timeF(l)=NaN;
                end
         end
        
         % Extract raster data set
         for l=1:length(timeF)
             if isnan(timeF(l))==0
                data(:,:,l)=eval(strcat('raster{',num2str(l),'}.',varF,'(:,:,1,',num2str(k),')'));
             else
                data(:,:,l)=NaN;
             end
         end
     
        
        % A. determine regression coefficients at cell levels
        % A.1. Identify grid locs for station data (both iMoMo Met Stations and WMO stations)
       
        
        
        gridLoc.WMO1 =  data(:,:,1);
        gridLoc.iMoMo1 =  data(:,:,1);

        baseV = 10^6;

        nWMO = 0;
        for idx = 1 : length(statData.WMO.LONLAT(:,1))
                temp = imbedm(statData.WMO.LONLAT(idx,2), statData.WMO.LONLAT(idx,1), ...
                    baseV + idx, gridLoc.WMO1, rasterR);
                gridLoc.WMO(:,:,idx) = temp;
                nWMO = nWMO + 1;
        end

        niMoMo = 0;
        for idx = 1 : length(statData.iMoMo.LONLAT(:,1))
                temp = imbedm(statData.iMoMo.LONLAT(idx,2), statData.iMoMo.LONLAT(idx,1), ...
                    baseV + idx, gridLoc.iMoMo1, rasterR);
                gridLoc.iMoMo(:,:,idx) = temp;
                niMoMo = niMoMo + 1;
        end

        gridLoc.WMO = gridLoc.WMO - baseV; gridLoc.WMO(gridLoc.WMO<0) = NaN;
        gridLoc.iMoMo = gridLoc.iMoMo - baseV; gridLoc.iMoMo(gridLoc.iMoMo<0) = NaN;


        % A.2.1 preparing time comparison (delete NaN values in WMO / iMoMo station data)
        for idx = 1 : length(statData.WMO.ID)
                statData.WMO.timeStatData{idx} = datenum(num2str(statData.WMO.YEARMODA(~isnan(eval(strcat('statData.WMO.',var2BC{1},'(:,idx)'))),idx)),'yyyymmdd');
                statData.WMO.idxNNan{idx} = ~isnan(eval(strcat('statData.WMO.',var2BC{1},'(:,idx)')));
                eval(['statData.WMO.TDataOnly{' num2str(idx) '} = statData.WMO.' var2BC{1} '(statData.WMO.idxNNan{' num2str(idx) '},' num2str(idx) ');'])
                %statData.WMO.TDataOnly{idx} = statData.WMO.TEMP(statData.WMO.idxNNan{idx},idx);
        end
        
        for idx = 1 : length(statData.iMoMo.ID)
                statData.iMoMo.timeStatData{idx} = datenum(num2str(statData.iMoMo.YEARMODA(~isnan(eval(strcat('statData.iMoMo.',var2BC{1},'(:,idx)'))),idx)),'yyyymmdd');
                statData.iMoMo.idxNNan{idx} = ~isnan(eval(strcat('statData.iMoMo.',var2BC{1},'(:,idx)')));
                eval(['statData.iMoMo.TDataOnly{' num2str(idx) '} = statData.iMoMo.' var2BC{1} '(statData.iMoMo.idxNNan{' num2str(idx) '},' num2str(idx) ');'])
                %statData.iMoMo.TDataOnly{idx} = statData.iMoMo.TEMP(statData.iMoMo.idxNNan{idx},idx);
        end
   
        % A.2.2 delete NAN entries in raster data
        for s = 1 : length(statData.WMO.ID)
            [idxGridR.WMO{s},idxGridC.WMO{s}] = find(gridLoc.WMO(:,:,s) == s);
            timeTNNaN.WMO{s} = timeF(~isnan(squeeze(data(idxGridR.WMO{s},idxGridC.WMO{s},:))));
            rasterNNaN.WMO{s} = squeeze(data(idxGridR.WMO{s},idxGridC.WMO{s},~isnan(squeeze(data(idxGridR.WMO{s},idxGridC.WMO{s},:)))));
        end

        for s = 1 : length(statData.iMoMo.ID)
            [idxGridR.iMoMo{s},idxGridC.iMoMo{s}] = find(gridLoc.iMoMo(:,:,s) == s);
            timeTNNaN.iMoMo{s} = timeF(~isnan(squeeze(data(idxGridR.iMoMo{s},idxGridC.iMoMo{s},:))));
            rasterNNaN.iMoMo{s} = squeeze(data(idxGridR.iMoMo{s},idxGridC.iMoMo{s},~isnan(squeeze(data(idxGridR.iMoMo{s},idxGridC.iMoMo{s},:)))));
        end
        
               
        % A.2.2 Matching data and times from stations and the grid data and getting
        % regression coeffs
        
      
        if debugging
            figure(2) % WMO
        end
        for s = 1 : length(statData.WMO.ID)
            [I,ia,ib] = intersect(statData.WMO.timeStatData{s},timeTNNaN.WMO{s});
            % regress
            A = statData.WMO.TDataOnly{s}(ia);
            B = rasterNNaN.WMO{s}(ib);
                      
            slope(k).WMO(s,1) = A\B; % then B = slope * A
            cC= corrcoef(A,B);     
            if  cC(1,2)<0.6
                slope(k).WMO(s,1)=1;
            end
            if slope(k).WMO(s,1)==0;
                slope(k).WMO(s,1)=1;
            end
                
            %Catch warning
            [warnmsg, msgid] = lastwarn;
            if strcmp(msgid,'MATLAB:rankDeficientMatrix')
               slope(k).WMO(s,1)=1;
               slope(k).error(s,1)=1;
            end
            warning('')
            
            
            if debugging
                % plotting
                subplot(3,1,s)
                scatter(statData.WMO.TDataOnly{s}(ia),rasterNNaN.WMO{s}(ib),'+')
                hold on
                plot([0 40],[0 40 * slope(k).WMO(s)],'k')
                hold off
                xlabel('WMO station data (T)'); ylabel('Gridded Data (GDAS T)')
                axis([0 40 0 40])
                grid minor
            end
        end

        if debugging
            figure(3) % iMoMo
        end
        
            
        
        for s = 1 : length(statData.WMO.ID)
            [I,ia,ib] = intersect(statData.iMoMo.timeStatData{s},timeTNNaN.iMoMo{s});
            % regress
            A = statData.iMoMo.TDataOnly{s}(ia);
            B = rasterNNaN.iMoMo{s}(ib);
            
         
            slope(k).iMoMo(s,1) = A\B; % then B = slope * A
            cC= corrcoef(A,B);     
            if  cC(1,2)<0.6
                slope(k).iMoMo(s,1)=1;
            end
            if slope(k).iMoMo(s,1)==0;
                slope(k).iMoMo(s,1)=1;
            end
            
           % Catch warning
            [warnmsg, msgid] = lastwarn;
            if strcmp(msgid,'MATLAB:rankDeficientMatrix')
               slope(k).iMoMo(s,1)=1;
               slope(k).error(s,1)=1;
            end
            warning('')
            
            if debugging
                % plotting
                subplot(3,1,s)
                scatter(statData.iMoMo.TDataOnly{s}(ia),rasterNNaN.iMoMo{s}(ib),'+')
                hold on
                plot([0 40],[0 40 * slope(k).iMoMo(s)],'k')
                hold off
                xlabel('iMoMo station data (T)'); ylabel('Gridded Data (GDAS T)')
                axis([0 40 0 40])
                grid minor
            end
        end
        
        % B. gIDW() for inverse distance weightning!
        % distance weights
        w = -1;
        % B.1. Define intrinsic grid coordinates
        XCoord = rasterR.XIntrinsicLimits(1):rasterR.XIntrinsicLimits(2)-1;
        YCoord = rasterR.YIntrinsicLimits(1):rasterR.YIntrinsicLimits(2)-1;
        [X,Y] = meshgrid(XCoord,YCoord); Xi = X(:); Yi = Y(:);
               
        
        iNNaNSlopeWMO = ~~slope(k).WMO; idxWMO.C = [idxGridC.WMO{:}]; idxWMO.R = [idxGridR.WMO{:}];
        iNNaNSlopeIMoMo = ~~slope(k).iMoMo; idxIMoMo.C = [idxGridC.iMoMo{:}]; idxIMoMo.R = [idxGridR.iMoMo{:}];
        Vi(k).WMO = gIDW(idxWMO.C(iNNaNSlopeWMO)'-0.5,idxWMO.R(iNNaNSlopeWMO)'-0.5,slope(k).WMO(iNNaNSlopeWMO),Xi,Yi,w);
        Vi(k).iMoMo = gIDW(idxIMoMo.C(iNNaNSlopeIMoMo)'-0.5,idxIMoMo.R(iNNaNSlopeIMoMo)'-0.5,slope(k).iMoMo(iNNaNSlopeIMoMo),Xi,Yi,w);

        Vi(k).WMO = reshape(Vi(k).WMO,size(data,1),size(data,2));
        Vi(k).iMoMo =reshape(Vi(k).iMoMo,size(data,1),size(data,2));
        
        % C. Return Values
        % Return a wighted average field where the weights correspond to how
        % trustworthy one believes the data is.
        Vi(k).weighted = statData.w.WMO * Vi(k).WMO + statData.w.iMoMo * Vi(k).iMoMo;
       
        if debugging
            figure(4)
            subplot(3,1,1)
            imagesc(Vi.WMO), colorbar
            title('WMO')
            subplot(3,1,2)
            imagesc(Vi.iMoMo), colorbar
            title('iMoMo')
            subplot(3,1,3)
            imagesc(Vi.weighted),colorbar
            title('Weigthed average')
        end
        
                
        
        %Reshape results to dimensions of input
        for l=1:length(timeF)
            if isfield(raster{l},varF)==1
                for m=1:sD
                    rasterDay(:,:,m,1)=eval(strcat('raster{',num2str(l),'}.',varF,'(:,:,',num2str(m),',',num2str(k),')./Vi(',num2str(k),').weighted'));
                end
                if strcmp(varF, 'FT')
                     rasterC{l}.FT(:,:,:,k)= rasterDay;
                elseif strcmp(varF, 'FP')
                     rasterC{l}.FP(:,:,:,k)= rasterDay;
                elseif strcmp(varF, 'FmaxT')
                     rasterC{l}.FmaxT(:,:,:,k)= rasterDay;
                elseif strcmp(varF, 'FminT')
                     rasterC{l}.FminT(:,:,:,k)= rasterDay;
                end
            end
        end
        
   end


elseif or(strcmp(dataSource,'GDAS'),strcmp(dataSource,'FEWS'))
    
    % 0. PLOTTING
    if debugging
        figure(1)
        geoshow(squeeze(raster(:,:,1)),rasterR,'DisplayType','texturemap'), colorbar
        hold on
        plot(statData.WMO.LONLAT(:,1),statData.WMO.LONLAT(:,2),'k*')
        plot(statData.iMoMo.LONLAT(:,1),statData.iMoMo.LONLAT(:,2),'r*')
        hold off
        legend('','WMO stations','iMoMo Stations')
    end
    
    % A. determine regression coefficients at cell levels
    % A.1. Identify grid locs for station data (both iMoMo Met Stations and WMO stations)
    gridLoc.WMO = squeeze(raster(:,:,1));
    gridLoc.iMoMo = squeeze(raster(:,:,1));
    
    baseV = 10^6;
    
    nWMO = 0;
    for idx = 1 : length(statData.WMO.LONLAT(:,1))
        temp = imbedm(statData.WMO.LONLAT(idx,2), statData.WMO.LONLAT(idx,1), ...
            baseV + idx, gridLoc.WMO, rasterR);
        gridLoc.WMO = temp;
        nWMO = nWMO + 1;
    end
    
    
    niMoMo = 0;
    for idx = 1 : length(statData.iMoMo.LONLAT(:,1))
        temp = imbedm(statData.iMoMo.LONLAT(idx,2), statData.iMoMo.LONLAT(idx,1), ...
            baseV + idx, gridLoc.iMoMo, rasterR);
        gridLoc.iMoMo = temp;
        niMoMo = niMoMo + 1;
    end
    
    gridLoc.WMO = gridLoc.WMO - baseV; gridLoc.WMO(gridLoc.WMO<0) = NaN;
    gridLoc.iMoMo = gridLoc.iMoMo - baseV; gridLoc.iMoMo(gridLoc.iMoMo<0) = NaN;
    
    if debugging
        keyboard
    end
    
    % A.2.1 preparing time comparison (delete NaN values in WMO / iMoMo station data)
    for idx = 1 : length(statData.WMO.ID)
        statData.WMO.timeStatData{idx} = datenum(num2str(statData.WMO.YEARMODA(~isnan(eval(strcat('statData.WMO.',var2BC{1},'(:,idx)'))),idx)),'yyyymmdd');
        statData.WMO.idxNNan{idx} = ~isnan(eval(strcat('statData.WMO.',var2BC{1},'(:,idx)')));
        eval(['statData.WMO.TDataOnly{' num2str(idx) '} = statData.WMO.' var2BC{1} '(statData.WMO.idxNNan{' num2str(idx) '},' num2str(idx) ');'])
        %statData.WMO.TDataOnly{idx} = statData.WMO.TEMP(statData.WMO.idxNNan{idx},idx);
    end
        
    for idx = 1 : length(statData.iMoMo.ID)
        statData.iMoMo.timeStatData{idx} = datenum(num2str(statData.iMoMo.YEARMODA(~isnan(eval(strcat('statData.iMoMo.',var2BC{1},'(:,idx)'))),idx)),'yyyymmdd');
        statData.iMoMo.idxNNan{idx} = ~isnan(eval(strcat('statData.iMoMo.',var2BC{1},'(:,idx)')));
        eval(['statData.iMoMo.TDataOnly{' num2str(idx) '} = statData.iMoMo.' var2BC{1} '(statData.iMoMo.idxNNan{' num2str(idx) '},' num2str(idx) ');'])
        %statData.iMoMo.TDataOnly{idx} = statData.iMoMo.TEMP(statData.iMoMo.idxNNan{idx},idx);
    end
    
    % A.2.2 delete NAN entries in raster data
    for s = 1 : length(statData.WMO.ID)
        [idxGridR.WMO{s},idxGridC.WMO{s}] = find(gridLoc.WMO == s);
        timeTNNaN.WMO{s} = timeT(~isnan(squeeze(squeeze(raster(idxGridR.WMO{s},idxGridC.WMO{s},:)))));
        rasterNNaN.WMO{s} = ...
            squeeze(raster(idxGridR.WMO{s},idxGridC.WMO{s},...
            ~isnan(squeeze(squeeze(raster(idxGridR.WMO{s},idxGridC.WMO{s},:))))));
    end
    
    for s = 1 : length(statData.iMoMo.ID)
        [idxGridR.iMoMo{s},idxGridC.iMoMo{s}] = find(gridLoc.iMoMo == s);
        timeTNNaN.iMoMo{s} = timeT(~isnan(squeeze(squeeze(raster(idxGridR.iMoMo{s},idxGridC.iMoMo{s},:)))));
        rasterNNaN.iMoMo{s} = ...
            squeeze(raster(idxGridR.iMoMo{s},idxGridC.iMoMo{s},...
            ~isnan(squeeze(squeeze(raster(idxGridR.iMoMo{s},idxGridC.iMoMo{s},:))))));
    end
    
    % A.2.2 Matching data and times from stations and the grid data and getting
    % regression coeffs
    if debugging
        figure(2) % WMO
    end
   
    
    for s = 1 : length(statData.WMO.ID)
        [I,ia,ib] = intersect(statData.WMO.timeStatData{s},timeTNNaN.WMO{s});
        % regress
        A = statData.WMO.TDataOnly{s}(ia);
        B = rasterNNaN.WMO{s}(ib);
        slope.WMO(s,1) = A\B; % then B = slope * A
        cC= corrcoef(A,B);     
            if  cC(1,2)<0.6
                slope(k).WMO(s,1)=1;
            end
            if slope.WMO(s,1)==0;
                slope.WMO(s,1)=1;
            end
        
        
        if debugging
            % plotting
            subplot(3,1,s)
            scatter(statData.WMO.TDataOnly{s}(ia),rasterNNaN.WMO{s}(ib),'+')
            hold on
            plot([0 40],[0 40 * slope.WMO(s)],'k')
            hold off
            xlabel('WMO station data (T)'); ylabel('Gridded Data (GDAS T)')
            axis([0 40 0 40])
            grid minor
        end
    end
    
       
    if debugging
        figure(3) % iMoMo
    end
    for s = 1 : length(statData.WMO.ID)
        [I,ia,ib] = intersect(statData.iMoMo.timeStatData{s},timeTNNaN.iMoMo{s});
        % regress
        A = statData.iMoMo.TDataOnly{s}(ia);
        B = rasterNNaN.iMoMo{s}(ib);
        slope.iMoMo(s,1) = A\B; % then B = slope * A
        
        cC= corrcoef(A,B);     
            if  cC(1,2)<0.6
                slope(k).iMoMo(s,1)=1;
            end
            if slope.iMoMo(s,1)==0;
                slope.iMoMo(s,1)=1;
            end
        
        
        if debugging
            % plotting
            subplot(3,1,s)
            scatter(statData.iMoMo.TDataOnly{s}(ia),rasterNNaN.iMoMo{s}(ib),'+')
            hold on
            plot([0 40],[0 40 * slope.iMoMo(s)],'k')
            hold off
            xlabel('iMoMo station data (T)'); ylabel('Gridded Data (GDAS T)')
            axis([0 40 0 40])
            grid minor
        end
    end
    
        
    % B. gIDW() for inverse distance weightning!
    % distance weights
    w = -1;
    % B.1. Define intrinsic grid coordinates
    XCoord = rasterR.XIntrinsicLimits(1):rasterR.XIntrinsicLimits(2)-1;
    YCoord = rasterR.YIntrinsicLimits(1):rasterR.YIntrinsicLimits(2)-1;
    [X,Y] = meshgrid(XCoord,YCoord); Xi = X(:); Yi = Y(:);
    iNNaNSlopeWMO = ~~slope.WMO; idxWMO.C = [idxGridC.WMO{:}]; idxWMO.R = [idxGridR.WMO{:}];
    iNNaNSlopeIMoMo = ~~slope.iMoMo; idxIMoMo.C = [idxGridC.iMoMo{:}]; idxIMoMo.R = [idxGridR.iMoMo{:}];
    Vi.WMO = gIDW(idxWMO.C(iNNaNSlopeWMO)'-0.5,idxWMO.R(iNNaNSlopeWMO)'-0.5,slope.WMO(iNNaNSlopeWMO),Xi,Yi,w);
    Vi.iMoMo = gIDW(idxIMoMo.C(iNNaNSlopeIMoMo)'-0.5,idxIMoMo.R(iNNaNSlopeIMoMo)'-0.5,slope.iMoMo(iNNaNSlopeIMoMo),Xi,Yi,w);
    
    Vi.WMO = reshape(Vi.WMO,size(raster,1),size(raster,2));
    Vi.iMoMo =reshape(Vi.iMoMo,size(raster,1),size(raster,2));
    
    % C. Return Values
    % Return a wighted average field where the weights correspond to how
    % trustworthy one believes the data is.
    
    Vi.weighted = statData.w.WMO * Vi.WMO + statData.w.iMoMo * Vi.iMoMo;
    rasterC = raster ./ repmat(Vi.weighted, [1 1 size(raster,3)]); % bias correction done!
    
    if debugging
        figure(4)
        subplot(3,1,1)
        imagesc(Vi.WMO), colorbar
        title('WMO')
        subplot(3,1,2)
        imagesc(Vi.iMoMo), colorbar
        title('iMoMo')
        subplot(3,1,3)
        imagesc(Vi.weighted),colorbar
        title('Weigthed average')
    end
    
    if debugging
        keyboard
    end
    
    
elseif strcmp(dataSource,'FEWS')
    
    
    
end

end % function