function [output,message] = testReadingSubAndForecastDataAndWritingDB()
%% function [output,message] = testReadingSubAndForecastDataAndWritingDB()
% 
% Reads in data stored under /data/processed/sub/ and data stored under
% results/forecast/ and compares them to the last file written to
% /data/processed/db. 
%
% @return output (boolean) 1 for success and 0 for failure.
% @return message (string) empty in case of success and error message in
%                          case of failure. 

output = 1;  % Success.
message = '';

clc
fprintf('----\n');
fprintf('Starting test reading sub and forecast data and writing of db.\n');
fprintf('----\n');

mfilepath = mfilename('fullpath');
mfilepathstr = strsplit(mfilepath,filesep);
current_dir = strjoin(mfilepathstr(1:end-1),filesep);
model_dir = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,'Themi'];

% Load latest ddddd.mat file.
dataBase_dir = [model_dir,filesep,'data',filesep,'processed'];
DB = load_data(dataBase_dir,'db');

% Load forecast data from results.
forecast_dir = [model_dir,filesep,'results',filesep,'forecast'];
FET = load_data(forecast_dir,'ET');
FG = load_data(forecast_dir,'G');
FQ = load_data(forecast_dir,'Q');
FS = load_data(forecast_dir,'S');

% Load forecast data from `sub`.
sub_dir = [model_dir,filesep,'data',filesep,'processed',filesep,'sub'];
FP = load_sub_data(sub_dir,'sub_FP.mat');
FT = load_sub_data(sub_dir,'sub_FT.mat');

% keyboard

% Tests.
bol = test_sub_data(DB,FP.FPsub{end},'P');
if bol ~= 1; output = 0; message = [message,'Failed in test for P.']; end
bol = test_sub_data(DB,FT.FTsub{end},'T');
if bol ~= 1; output = 0; message = [message,'Failed in test for T.']; end

bol = test_data(DB,FET.FET.catch,'ET');
if bol ~= 1; output = 0; message = [message,'Failed in test for ET.']; end
bol = test_data(DB,FG.FG.catch,'G');
if bol ~= 1; output = 0; message = [message,'Failed in test for G.']; end
bol = test_data(DB,FQ.FQ.catch,'Q');
if bol ~= 1; output = 0; message = [message,'Failed in test for Q.']; end
bol = test_data(DB,FS.FS.catch,'S');
if bol ~= 1; output = 0; message = [message,'Failed in test for S.']; end

fprintf('DONE.\n');

end


%% Local functions

function data = load_data(directory,name)
dir_list = dir([directory,filesep,name]);
index = length(dir_list);
file = cell2mat({dir_list(end).name}');  % Get the latest file.
while (index > 1 && (isempty(strfind(file,'.mat')) || ...
                     ~isempty(regexp(file,'^\.','once')) || ...
                     ~isempty(regexp(file,'^\d{6}_\w{2,3}_orig.mat$','once'))))
  index = index - 1;
  file = cell2mat({dir_list(index).name}');
end
try
  data = load([directory,filesep,name,filesep,file]);
  fprintf('Loading %s.\n',file);
catch
  data = [];
  fprintf('   ERROR: Problem loading %s from %s.\n',name,[directory,filesep,name,filesep,file])
end
end


function data = load_sub_data(directory,filename)
try
  data = load([directory,filesep,filename]);
  fprintf('Loading %s.\n',filename);
catch
  data = [];
  fprintf('   ERROR: Problem loading %s from %s.\n',filename,directory)
end
end



function out = test_sub_data(data_base,forecast,name)
out = 1; % Success;
% Compute stats.
fname = ['F',name];
if isfield(forecast,['F','min',name'])
  data_min = min(forecast.(['Fmin',name]),[],3);
  data_max = max(forecast.(['Fmax',name]),[],3);
else
  data_max = max(forecast.(fname),[],3);
  data_min = min(forecast.(fname),[],3);
end
data_mean = mean(forecast.(fname),3);
data_time = forecast.timeF;
date=cellstr(strcat(datestr(data_time,'yyyy'),'-',datestr(data_time,'mm'),'-',datestr(data_time,'dd'),' 12:00:00'));

% Test date
fprintf('Testing %s date . . . ',name);
if ~isequal(data_base.dateMet,date)
  fprintf('Failure at %s date.\n',name)
  out = 0;
else
  fprintf('. . . done.\n');
end

% Test mean, max, min
fprintf('Testing %s . . . ',name);
if ~isequal(data_base.(name)(:,:,1),data_mean)
  fprintf('Failure at %smean.\n',name);
  out = 0;
%   keyboard
end
if ~isequal(data_base.(name)(:,:,2),data_max)
  fprintf('Failure at %smax.\n',name)
  out = 0;
  keyboard
end
if ~isequal(data_base.(name)(:,:,3),data_min)
  fprintf('Failaure at %smin.\n',name)
  out = 0;
else
  fprintf('. . . done.\n');
end
end


function out = test_data(data_base,forecast,name)
out = 1;
nSubCatchments = size(forecast,2);
fprintf('Testing %s . . . ',name);
for j=1:nSubCatchments
  for k=1:length(forecast{j}.Q) % Ensemble size
    maxMQ(:,k)=max(forecast{j}.Q{k},[],2);
    minMQ(:,k)=min(forecast{j}.Q{k},[],2);
    meanMQ(:,k)=mean(forecast{j}.Q{k},2);
  end
Q(:,j,1)=mean(meanMQ,2);
Q(:,j,2)=max(maxMQ,[],2);
Q(:,j,3)=min(minMQ,[],2);
end
if ~isequal(data_base.(name),Q)
  fprintf('Failure at %s.\n',name)
  out = 0;
else
  fprintf('. . . done.\n');
end
  
end