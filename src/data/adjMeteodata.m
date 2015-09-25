function [cData] = adjMeteodata(data, param, option, DEM, iC, iE)

%% ADJUST METEOROLOGICAL DATA
%Calculate Correction Factor


switch option
    case 1 %Polynominal factor correction
        cF=1+param{1}.k1*(DEM(iE)./1000-DEM(iC)./1000)+param{1}.k2*(DEM(iE)./1000-DEM(iC)./1000).^2;
        % Adjust cells above thresholds
        for l=1:size(data,3)
            lP=data(:,:,l);
            lP(iE)=lP(iC).*cF;
            cData(:,:,l)=reshape(lP,size(data,1),size(data,2));
        end
           
    case 2 %Additive lapse rate
        cF=(DEM(iE)-DEM(iC)).*param{2}.lapse/100;
        % Adjust cells above thresholds
        for l=1:size(data,3)
            lT=data(:,:,l);
            lT(iE)=lT(iC)+cF;
            cData(:,:,l)=reshape(lT,size(data,1),size(data,2));
        end
end
  

end