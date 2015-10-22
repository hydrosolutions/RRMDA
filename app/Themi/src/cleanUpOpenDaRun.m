function [] = cleanUpOpenDaRun()
%% function [] = cleanUpOpenDaRun()
%
% Load all odaE.mat and odaS0G0.mat files from the model/output/workX
% directories and concatenate the data.

fid = fopen('cleanUpOpenDaRun.log','w+');
doublefprintf(fid,'%s Starting cleanUpOpenDaRun.\n',datestr(now));

try
  
  doublefprintf(fid,'Setting up.\n');
  mfilepath = mfilename('fullpath');
  mfilepathcell = strsplit(mfilepath,filesep);
  currentDir = strjoin(mfilepathcell(1:end-1),filesep); % Remove file name.
  outputDir = [currentDir,filesep,'model',filesep,'output'];
  
  %% Load oda files from model/output/workX directories.
  doublefprintf(fid,'Loading work files.\n');
  try
  list = dir(outputDir);
  numberOfDirsInList = length(list);
  E = [];
  S0 = [];
  G0 = [];
  SG = [];
  index = 1;
  for i = 1:numberOfDirsInList  %// Iterate over all entries in outputDir.
    if ~isempty(strfind(list(i).name,'work'))  %// The directory is a working directory.
      try
        load([outputDir,filesep,list(i).name,filesep,'odaE.mat']);
      catch le
        doublefprintf(fid,'Problem loading %s.\n',[outputDir,filesep,flist(i).name,filesep,'odaE.mat']);
        rethrow(le)
      end
      odaE(index).E = [];
      odaE(index).E = E;
      try
        load([outputDir,filesep,list(i).name,filesep,'odaS0G0.mat']);
      catch le
        doublefprintf(fid,'Problem loading %s.\n',[outputDir,filesep,flist(i).name,filesep,'odaS0G0.mat']);
        rethrow(le)
      end
      odaS(index).S0 = S0;
      odaS(index).G0 = G0;
      odaS(index).SG = SG;
      index = index + 1;
    end
  end
  catch e
    doublefprintf(fid,'Problem loading work files.\n');
    rethrow(e)
  end
  
  clear E S0 G0 SG
  
  %% Concatenate the content of the files.
  doublefprintf(fid,'Concatenating data.\n');
  try
  E = [];
  S0 = odaS(1).S0;
  G0 = odaS(1).G0;
  SG = [];
  %// Take first row of work data and store it in output variable.
  numberOfWorkingDirs = length(odaE);
  for i = 1:numberOfWorkingDirs
    E(:,i) = odaE(i).E(:,1);  
    SG(:,i) = odaS(i).SG(:,1);
  end
  catch e
    doublefprintf(fid,'Problem concatenating data.\n');
    rethrow(e)
  end
  
  % Read in original input files.
  origE = load([currentDir,filesep,'model',filesep,'input',filesep,'E.mat']);
  origS = load([currentDir,filesep,'model',filesep,'input',filesep,'S0G0.mat']);
  
  for i = size(E,2)+1:size(origE.E,2)
    E(:,i) = origE.E(:,i);
    SG(:,i) = origS.SG(:,i);
  end
  
  %% Saving data.
  doublefprintf(fid,'Saving the data.\n');
  try
  save([currentDir,filesep,'model',filesep,'input',filesep,'E.mat'],'E');
  save([currentDir,filesep,'model',filesep,'input',filesep,'S0G0.mat'],'S0','G0','SG');
  catch e
    doublefprintf(fid,'Problem saving data.\n');
    rethrow(e)
  end
  
  doublefprintf(fid,'%s DONE.\n',datestr(now));
  fclose(fid);
  %exit(0)

catch me
  
  %doublefprintf(fid,'Problem in cleanUpOpenDaRun.\n');
  fclose(fid);
  rethrow(me)
  %exit(1)

end

end


function doublefprintf(fid,message,varargin)
%% function doubpefprintf prints to file handle and to console.
%
% function doublefprintf(fid,message,varargin)
%
% @param fid (int) file handle to a writable file.
% @param message (string) with message to print.
% @param varargin (any) optional arguments to be passed to the message.
%

if nargin < 3
    fprintf(fid,message);
    fprintf(message);
elseif nargin >= 3
    fprintf(fid,message,varargin{:});
    fprintf(message,varargin{:});
end

end
