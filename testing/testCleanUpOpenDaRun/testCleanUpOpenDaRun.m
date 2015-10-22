function [] = testCleanUpOpenDaRun()
%% function [] = testCleanUpOpenDaRun()
%
% Tests the function cleanUpOpenDaRun.m in a test environment separate from
% runModel. 
%


ensembleSize = 10;

warning('off','MATLAB:rmpath:DirNotFound');  %// Turn off unneccessary warning.

%% Clean up if neccessary.
if exist(['testEnvironment',filesep,'model'],'dir') == 7
  rmdir(['testEnvironment',filesep,'model'],'s'); % Remove directory and sub-directories.
end

if exist(['testEnvironment',filesep,'cleanUpOpenDaRun.m'],'file') == 2
  delete(['testEnvironment',filesep,'cleanUpOpenDaRun.m']);
end

if exist(['testEnvironment',filesep,'model',filesep,'input'],'dir') == 7
  rmdir(['testEnvironment',filesep,'model',filesep,'input'],'s');
end

%% Set up test environment.
try

  % Get name of current directory.
  mfilepath = mfilename('fullpath');
  mfilepathcell = strsplit(mfilepath,filesep);
  current_dir = strjoin(mfilepathcell(1:end-1),filesep); % Remove file name.
  
  try
    fid = fopen(strcat(current_dir,filesep,'testCleanUpOpenDaRun.log'),'w+');
    doublefprintf(fid,'-----\n');
    c = clock;
    doublefprintf(fid,'START testCleanUpOpenDaRun() %2d.%2d %4d, %2d:%2d\n',c(3),c(2),c(1),c(4),c(5));
  catch error
    fprintf('Problem opening file log for writing.\n')
    fprintf('Message: %s\n',error.message);
    rethrow(error)
  end
  doublefprintf(fid,'Setting up test environment. . . ');
  
  % Copy most recent cleanUpOpenDaRun.m from RRMDA/src/setup.
  src_file = [current_dir,filesep,'..',filesep,'..',filesep,'src',filesep,'setup',filesep,'cleanUpOpenDaRun.m'];
  dest_dir = [current_dir,filesep,'testEnvironment',filesep];
  copyfile(src_file,dest_dir);
  
  % Add testEnvironment to MATLABPATH.
  %addpath(genpath([current_dir,filesep,'testEnvironment']));
  
  % Set up workX directories.
  mkdir([dest_dir,'model',filesep,'output']);
  mkdir([dest_dir,'model',filesep,'input']);
  copyfile([dest_dir,'testData',filesep,'E.mat'],[dest_dir,'model',filesep,'input',filesep]);
  copyfile([dest_dir,'testData',filesep,'S0G0.mat'],[dest_dir,'model',filesep,'input',filesep]);
  odaE = load([dest_dir,'testData',filesep,'E.mat']);
  odaS0G0 = load([dest_dir,'testData',filesep,'S0G0.mat']);
  for i=1:ensembleSize
    workingDir = [dest_dir,'model',filesep,'output',filesep,strcat('work',num2str(i-1))]; 
    mkdir(workingDir)
    E = odaE.E(:,i);
    S0G0.SG = odaS0G0.SG(:,i);
    S0G0.S0 = odaS0G0.S0;
    S0G0.G0 = odaS0G0.G0;
    save([workingDir,filesep,'odaE.mat'],'E');
    save([workingDir,filesep,'odaS0G0.mat'],'-struct','S0G0');
  end
  
  doublefprintf(fid,'. . . done.\n');

  %% Call to cleanUpOpenDaRun
  doublefprintf(fid,'Calling cleanUpOpenDaRun. . . ');
  cd(dest_dir);
  cleanUpOpenDaRun()
  cd(current_dir);
  doublefprintf(fid,'. . . done.\n');
  
  %% Compare result to expected data loaded above (odaE and odaS0G0).
  doublefprintf(fid,'Testing resutls . . . ');
  E = [];
  S0 = [];
  G0 = [];
  SG = [];
  load([dest_dir,'model',filesep,'input',filesep,'E.mat']);
  load([dest_dir,'model',filesep,'input',filesep,'S0G0.mat']);
  
  if ~isequal(E,odaE.E)
    doublefprintf(fid,'Test E failed.\n');
  end
  
  if ~isequal(S0,odaS0G0.S0)
    doublefprintf(fid,'Test S0 failed.\n');
  end
  
  if ~isequal(G0,odaS0G0.G0)
    doublefprintf(fid,'Test G0 failed.\n');
  end
  
  if ~isequal(SG,odaS0G0.SG)
    doublefprintf(fid,'Test SG failed.\n');
  end

  doublefprintf(fid,'. . . done.\n');
  
  %% Clean up again.
  doublefprintf(fid,'Cleaning up test run . . . ');
  %rmpath(genpath([current_dir,filesep,'testEnvironment']));
  
  if exist(['testEnvironment',filesep,'model'],'dir') == 7
    rmdir(['testEnvironment',filesep,'model'],'s'); % Remove directory and sub-directories.
  end

  if exist(['testEnvironment',filesep,'cleanUpOpenDaRun.m'],'file') == 2
    delete(['testEnvironment',filesep,'cleanUpOpenDaRun.m']);
  end
  
  if exist(['testEnvironment',filesep,'model',filesep,'input'],'dir') == 7
    rmdir(['testEnvironment',filesep,'model',filesep,'input'],'s');
  end
  
  doublefprintf(fid,'. . . done.\n');
  doublefprintf(fid,'DONE testing.\n');
  
catch e
  
  cd(current_dir)
  fclose('all');
  rethrow(e)
  
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