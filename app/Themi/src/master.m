% Operational script to execute Themi RRM Model

% VERSION: 24/04/2015

% Copyright (C) 2015 hydrosolutions ltd. Zurich, Switzerland

gR=timer('TimerFcn','getRaw_Themi');        
pR=timer('TimerFcn','processRaw_Themi');
rM=timer('TimerFcn','runModel_Themi');
sD=timer('TimerFcn','sendtoDB_Themi');
mM=timer('TimerFcn','MatlabMail_Themi');
cD=timer('TimerFcn','controlData_Themi');

set(gR, 'ExecutionMode', 'fixedRate') 
set(gR, 'Period', 10800); 
set(pR, 'ExecutionMode', 'fixedRate') 
set(pR, 'Period', 10800); 
set(rM, 'ExecutionMode', 'fixedRate') 
set(rM, 'Period', 10800); 
set(sD, 'ExecutionMode', 'fixedRate') 
set(sD, 'Period', 10800);
set(mM, 'ExecutionMode', 'fixedRate') 
set(mM, 'Period', 86400); 
set(cD, 'ExecutionMode', 'fixedRate') 
set(cD, 'Period', 86400); 

% start times crucial to have data from previous day available!
startat(gR, '17:30:00'); 
startat(pR, '17:40:00'); 
startat(rM, '17:50:00'); 
startat(sD, '18:00:00');
startat(mM, '18:15:00'); 
startat(cD, '18:20:00'); 



