function [out,message] = testAssimilation()
%% function [out,message] = testAssimilation()
%
% Tests the function assimilation.m in a test environment separate from
% runModel. 
%

out = 1;
message = '';


%% Clean up if neccessary.
if exist('testEnvironment/src','dir') == 7
  rmdir('testEnvironment/src','s'); % Remove directory and sub-directories.
end

%% Load test environment.
try
  load('testSetup.mat');
  
  % Get name of current directory.
  mfilepath = mfilename('fullpath');
  mfilepathcell = strsplit(mfilepath,filesep);
  current_dir = strjoin(mfilepathcell(1:end-1),filesep); % Remove file name.
  
  try
    fid = fopen(strcat(current_dir,filesep,'testAssimilation.log'),'w+');
    doublefprintf(fid,'-----\n');
    c = clock;
    doublefprintf(fid,'START testAssimilation() %2d.%2d %4d, %2d:%2d\n',c(3),c(2),c(1),c(4),c(5));
  catch error
    fprintf('Problem opening file log for writing.\n')
    fprintf('Message: %s\n',error.message);
    rethrow(error)
  end
  doublefprintf(fid,'Setting up test environment. . . ');
  
  % Copy RRM/src directory containing all data assimilation source codes.
  src_dir = [current_dir,filesep,'..',filesep,'..',filesep,'src'];
  dest_dir = [current_dir,filesep,'testEnvironment',filesep,'src',filesep];
  copyfile(src_dir,dest_dir);
  
  sPath = setup.Path(1:end-4);
  
  % Copy original restart files.
  restart_dir = ['testEnvironment',filesep,'app',filesep,'TestThemi',filesep,...
                 'resources',filesep,'restart'];
  copyfile([restart_dir,filesep,'test_E.mat'],[restart_dir,filesep,'E.mat']);
  copyfile([restart_dir,filesep,'test_S0G0.mat'],[restart_dir,filesep,'S0G0.mat']);
  
  % Add paths of test environment to MATLABPATH.
  addpath(genpath([sPath,filesep,'src']));
  addpath(genpath(setup.mPath));
  
  paths.home = strcat(sPath,'/src/enkf/');
  paths.main = strcat(setup.mPath,'/');
  paths.data = strcat(setup.mPath,'/data/');
  
  cd(paths.home)
  doublefprintf(fid,'. . . done.\n');
catch e
  message = [message,sprintf('Problem setting up test environment. Error message: %s. Returning.\n',e.message)];
  out = 0;
  return
end
%% Call assimilation.
doublefprintf(fid,'Calling assimilation.m . . . ');
try
  assimilation(paths);
  doublefprintf(fid,'. . . done.\n')
catch e
  message = [message,sprintf('Problem in call to assimilation(). Error message: %s. Returning.\n',e.message)];
  doublefprintf(fid,'Problem: %s.\n',e.message);
  out = 0;
  return
end
  
cd(current_dir);

%% Test against expected values. 
result_dir = [current_dir,filesep,'testEnvironment',filesep,'app',filesep,...
              'TestThemi',filesep,'results'];

file = '736237_E.mat';
doublefprintf(fid,'Testing %s . . . ',file);
[bol,text] = testData([result_dir,filesep,'nowcast',filesep,file],...
                      ['expectedResults',filesep,file]);
if bol ~= 1
  message = [message,sprintf('Failed in test for %s.\nMessage: %s.\n',file,text)];
  out = bol;
  doublefprintf(fid,'Problem: %s.\n',message);
else
  doublefprintf(fid,'. . . done.\n');
end

file = '736239_FET.mat';
doublefprintf(fid,'Testing %s . . . ',file);
[bol,text] = testData([result_dir,filesep,'forecast',filesep,'ET',filesep,file],...
                      ['expectedResults',filesep,file]);
if bol ~= 1
  message = [message,sprintf('Failed in test for %s.\nMessage: %s.\n',file,text)];
  out = bol;
  doublefprintf(fid,'Problem: %s.\n',message);
else
  doublefprintf(fid,'. . . done.\n');
end

file = '736239_FG.mat';
doublefprintf(fid,'Testing %s . . . ',file);
[bol,text] = testData([result_dir,filesep,'forecast',filesep,'G',filesep,file],...
                      ['expectedResults',filesep,file]);
if bol ~= 1
  message = [message,sprintf('Failed in test for %s.\nMessage: %s.\n',file,text)];
  out = bol;
  doublefprintf(fid,'Problem: %s.\n',message);
else
  doublefprintf(fid,'. . . done.\n');
end

file = '736239_FQ.mat';
doublefprintf(fid,'Testing %s . . . ',file);
[bol,text] = testData([result_dir,filesep,'forecast',filesep,'Q',filesep,file],...
                      ['expectedResults',filesep,file]);
if bol ~= 1
  message = [message,sprintf('Failed in test for %s.\nMessage: %s.\n',file,text)];
  out = bol;
  doublefprintf(fid,'Problem: %s.\n',message);
else
  doublefprintf(fid,'. . . done.\n');
end

file = '736239_FS.mat';
doublefprintf(fid,'Testing %s . . . ',file);
[bol,text] = testData([result_dir,filesep,'forecast',filesep,'S',filesep,file],...
                      ['expectedResults',filesep,file]);
if bol ~= 1
  message = [message,sprintf('Failed in test for %s.\nMessage: %s.\n',file,text)];
  out = bol;
  doublefprintf(fid,'Problem: %s.\n',message);
else
  doublefprintf(fid,'. . . done.\n');
end


%% Clean up.
try
  rmpath(genpath([sPath,filesep,'src']));
  rmpath(genpath(setup.mPath));
  
  if exist('testEnvironment/src','dir') == 7
    rmdir('testEnvironment/src','s'); % Remove directory and sub-directories.
  end
  
  deleteIfExists([result_dir,filesep,'nowcast',filesep,'736237_E.mat']);
  deleteIfExists([result_dir,filesep,'forecast',filesep,'ET',filesep,'736239_FET.mat']);
  deleteIfExists([result_dir,filesep,'forecast',filesep,'G',filesep,'736239_FG.mat']);
  deleteIfExists([result_dir,filesep,'forecast',filesep,'Q',filesep,'736239_FQ.mat']);
  deleteIfExists([result_dir,filesep,'forecast',filesep,'S',filesep,'736239_FS.mat']);
  
  deleteIfExists(['testEnvironment',filesep,'app',filesep,'TestThemi',filesep,...
                  'resources',filesep,'restart',filesep,'E.mat']);
  deleteIfExists(['testEnvironment',filesep,'app',filesep,'TestThemi',filesep,...
                  'resources',filesep,'restart',filesep,'S0G0.mat']);
catch e
  message = [message,sprintf('Problem cleaning up after tests. Error message: %s. Returning.\n',e.message)];
  out = 0;
  return
end

doublefprintf(fid,'DONE.\n');
fclose(fid);

end


%% Local functions.

function [outbol,outtext] = testData(computedDataFileName,expectedDataFileName)

outbol = 1;
outtext = '';

try
  expectedData = load(expectedDataFileName); 
catch me
  outtext = sprintf('Problem loading %s. Error message: %s.\n',...
                    expectedDataFileName,me.message);
  outbol = 0;
  return
end
try
  computedData = load(computedDataFileName); 
catch me
  outtext = sprintf('Problem loading %s. Error message: %s.\n',...
                    expectedDataFileName,me.message);
  outbol = 0;
  return
end
if ~isequal(computedData,expectedData)
  outtext = 'Error in test: computed data ~= expected data.\n';
  outbol = 0;
  return
end

end


function [] = deleteIfExists(fileName)

if exist(fileName,'file') == 2
  delete(fileName);
end

end


function doublefprintf(fid,message,varargin)
%% function doubpefprintf prints to file handle and to console.
%
% function doublefprintf(fid,message,varargin)
%
% @param fid int file handle to a writable file.
% @param message string with message to print.
% @param varargin arguments to be passed to the message.
%

if nargin < 3
    fprintf(fid,message);
    fprintf(message);
elseif nargin >= 3
    fprintf(fid,message,varargin{:});
    fprintf(message,varargin{:});
end

end