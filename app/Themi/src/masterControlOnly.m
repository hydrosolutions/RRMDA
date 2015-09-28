% Operational script to execute Themi RRM Model

% VERSION: 24/04/2015

cD=timer('TimerFcn','controlData_Themi');


set(cD, 'ExecutionMode', 'fixedRate') 
set(cD, 'Period', 86400); 

% start times crucial to have data from previous day available!
 
startat(cD, '18:20:00'); 



