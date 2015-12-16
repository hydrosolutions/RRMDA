function [] = oda_adaptSetupMat()
% Is copied to the openDA working directory. There it gets address of the
% path to this file and adapts the paths in setup.mat accordingly. 

try
% Get current directory.
filepath = mfilename('fullpath');
[path,~,~] = fileparts(filepath);
parts = strsplit(path,filesep);
modelname = parts(end-1);

% Load original setup.mat
load('setup.mat');
delete('setup.mat');

l = length(parts)-2;
setup.mPath = '';
setup.Path = '';
for i = 2:l
  setup.mPath = [setup.mPath,filesep,char(parts(i))];
  setup.Path = [setup.Path,filesep,char(parts(i))];
end
setup.mPath = [setup.mPath,filesep,char(parts(l+1))];

save('setup.mat','setup');

exit(0)

catch ME
  
  fprintf('Problem in oda_adaptSetupMat()!\n')
  rethrow(ME)
  exit(1)



end