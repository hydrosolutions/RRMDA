% Operational script to execute Themi RRM Model
%
% AUTHOR: Tobias Siegfried, 02/06/2015, based on Sebastian Stoll script.
% 
% VERSION: 02/06/2015, added various timer options and changed the BusyMode
% option.

% Copyright (C) 2015 hydrosolutions ltd. Zurich, Switzerland

gR = timer('TimerFcn','getRaw_Themi');        
pR = timer('TimerFcn','processRaw_Themi');
rM = timer('TimerFcn','runModel_Themi');
sD = timer('TimerFcn','sendtoDB_Themi');
mM = timer('TimerFcn','MatlabMail_Themi');
cD = timer('TimerFcn','controlData_Themi');

set(gR,'ExecutionMode','fixedRate','BusyMode','queue','Name','getRaw','Period',10800) 
% gR.StartFcn = @(~,thisEvent)disp([thisEvent.Type ' getRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% gR.TimerFcn = @(~,thisEvent)disp([thisEvent.Type ' getRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% gR.StopFcn = @(~,thisEvent)disp([thisEvent.Type ' getRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);

set(pR,'ExecutionMode','fixedRate','BusyMode','queue','Name','process_Raw','Period',10800) 
% pR.StartFcn = @(~,thisEvent)disp([thisEvent.Type ' processRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% pR.TimerFcn = @(~,thisEvent)disp([thisEvent.Type ' processRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% pR.StopFcn = @(~,thisEvent)disp([thisEvent.Type ' processRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);

set(rM,'ExecutionMode','fixedRate','BusyMode','queue','Name','runModel','Period',10800) 
% rM.StartFcn = @(~,thisEvent)disp([thisEvent.Type ' processRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% rM.TimerFcn = @(~,thisEvent)disp([thisEvent.Type ' processRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% rM.StopFcn = @(~,thisEvent)disp([thisEvent.Type ' processRaw executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);

set(sD,'ExecutionMode','fixedRate','BusyMode','queue','Name','sendtoDB','Period',10800) 
% sD.StartFcn = @(~,thisEvent)disp([thisEvent.Type ' sendtoDB executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% sD.TimerFcn = @(~,thisEvent)disp([thisEvent.Type ' sendtoDB executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% sD.StopFcn = @(~,thisEvent)disp([thisEvent.Type ' sendtoDB executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);

set(mM,'ExecutionMode','fixedRate','BusyMode','queue','Name','MatlabMail','Period',86400) 
% mM.StartFcn = @(~,thisEvent)disp([thisEvent.Type ' matlabMail executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% mM.TimerFcn = @(~,thisEvent)disp([thisEvent.Type ' matlabMail executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% mM.StopFcn = @(~,thisEvent)disp([thisEvent.Type ' matlabMail executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);

set(cD,'ExecutionMode','fixedRate','BusyMode','queue','Name','controlData','Period',86400) 
% cD.StartFcn = @(~,thisEvent)disp([thisEvent.Type ' controlData executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% cD.TimerFcn = @(~,thisEvent)disp([thisEvent.Type ' controlData executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
% cD.StopFcn = @(~,thisEvent)disp([thisEvent.Type ' controlData executed ' datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);

% start times crucial to have data from previous day available!
startat(gR, '17:30:00'); 
startat(pR, '17:40:00'); 
startat(rM, '17:50:00'); 
startat(sD, '18:00:00');
startat(mM, '18:15:00'); 
startat(cD, '18:20:00'); 



