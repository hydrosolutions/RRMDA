function [out_args] = compare_data(in_args)
% FUNCTION [out_args] = compare_results(in_args)
%
% Reads data from data_for_comparison/ and compares them to results of the
% current offline model run if data files with the same name are available.
% Tests are done depending on the files present in data_for_comparison.
%
%
%
% Author Beatrice Marti, hydrosolutions ltd. (marti@hydrosolutions.ch)
%
% Copyright (C) 2015 hydrosolutions ltd. Zurich, Switzerland
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.

%% Global variables.
mfilepath = mfilename('fullpath');
mfilepathstr = strsplit(mfilepath,filesep);
current_dir = strjoin(mfilepathstr(1:end-1),filesep);
test_dir = [current_dir,filesep,'data_for_comparison'];

%% Open log file.
try
  fid = fopen(strcat(current_dir,filesep,'compare_data.log'),'w+');
  doublefprintf(fid,'-----\n');
  c = clock;
  doublefprintf(fid,'START compare_data() %2d.%2d %4d, %2d:%2d\n',c(3),c(2),c(1),c(4),c(5));
catch e
  fprintf('Problem opening file log for writing.\n')
  fprintf('Message: %s\n',e.message);
  rethrow(e)
end

%% Get file names.
% Get file names from data_for_comparison/
list_dir_data_for_comparison = dir(test_dir);

% If the directory is empty print a warning to the log file and return.
number_of_files = size(list_dir_data_for_comparison,1);
if number_of_files == 2
  doublefprintf(fid,'Warning: Directory %s is empty. No tests are done.\n',test_dir);
  return
end

%% Tests
for i = 3:number_of_files % Start at the third file since the first two are directories '.' and '..'.
  
  file = list_dir_data_for_comparison(i).name;
  doublefprintf(fid,'Testing file %s. . . ',file);
  
  % Ignore files starting with .
  if regexp(file,'^\.','once')
    doublefprintf(fid,'. . . file ignored.\n');
    continue
  end
  
  % Testing getRaw.
  if strcmp('FEWS.mat',file) % Returns 1 (true) or 0 (false)
    % Test if FEWS.mat is present in the local directory hierarchy.
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'raw',filesep,'FEWS.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(local_data.timeP,data.timeP)
        doublefprintf(fid,'\n   WARNING: timeP of %s not equal %s.\n',local_file,file);
      end
      if ~isequal(local_data.R_P,data.R_P)
        doublefprintf(fid,'\n   WARNING: R_P of %s not equal %s.\n',local_file,file);
      end
      if ~isequal(local_data.P,data.P)
        diff = local_data.P-data.P;
        if any(any(abs(diff) > 0.001))
          doublefprintf(fid,'\n   WARNING: P of %s nod equal %s within tolerance.\n',local_file,file);
        end
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
  
  elseif strcmp('GDAS.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'raw',filesep,'GDAS.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(local_data.timeT,data.timeT)
        doublefprintf(fid,'\n   WARNING: timeT of %s not equal %s.\n',local_file,file);
      end
      if ~isequal(local_data.R_T,data.R_T)
        doublefprintf(fid,'\n   WARNING: R_T of %s not equal %s.\n',local_file,file);
      end
      if ~isequal(local_data.T,data.T)
        diff = local_data.T-data.T;
        if any(any(abs(diff) > 0.001))
          doublefprintf(fid,'\n   WARNING: T of %s nod equal %s within tolerance.\n',local_file,file);
        end
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
  
  elseif strcmp('GFS.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'raw',filesep,'GFS.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(size(local_data.F),size(data.F))
        doublefprintf(fid,'\n   WARNING: size(F) of %s not equal %s.\n',local_file,file);
      end
      for j = size(local_data.F,2):-50:1
        if ~isequal(local_data.F{j}.timeF,data.F{j}.timeF)
          doublefprintf(fid,'\n   WARNING: F{%d}.timeF of %s not equal %s.\n',j,local_file,file);
        end
        if ~isequal(local_data.F{j}.FT,data.F{j}.FT)
          diff = local_data.F{j}.FT - data.F{j}.FT;
          if any(any(abs(diff) > 0.001))
            doublefprintf(fid,'\n   WARNING: F{%d}.FT of %s not equal to %s within tolerance.\n',j,local_file,file);
          end
        end
        if ~isequal(local_data.F{j}.FP,data.F{j}.FP)
          diff = local_data.F{j}.FP - data.F{j}.FP;
          if any(any(abs(diff) > 0.001))
            doublefprintf(fid,'\n   WARNING: F{%d}.FP of %s not equal %s within tolerance.\n',j,local_file,file);
          end
        end
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
  
  elseif strcmp('iMoMo.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'raw',filesep,'iMoMo.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(local_data.iMoMo.ID,data.iMoMo.ID)
        doublefprintf(fid,'\n   WARNING: ID of %s not equal %s.\n',local_file,file);
      end
      if ~isequal(size(local_data.iMoMo.TEMP,1),size(data.iMoMo.TEMP,1))
        doublefprintf(fid,'   ERROR in size(iMoMo.TEMP).\n');
        doublefprintf(fid,'         size(local_data.iMoMo.TEMP,1) = %d.\n',size(local_data.iMoMo.TEMP,1));
        doublefprintf(fid,'         size(      data.iMoMo.TEMP,1) = %d.\n',size(data.iMoMo.TEMP,1));
        return
      else
        if ~isequal(local_data.iMoMo.TEMP,data.iMoMo.TEMP)
          abs_diff = local_data.iMoMo.TEMP-data.iMoMo.TEMP;
          rel_diff1 = abs_diff ./local_data.iMoMo.TEMP;
          rel_diff2 = abs_diff ./data.iMoMo.TEMP;
          if any(any(abs(rel_diff1) > 0.05 & abs(rel_diff2) > 0.05))
            doublefprintf(fid,'\n   WARNING: TEMP of %s nod equal %s within tolerance.\n',local_file,file);
            keyboard
          end
        end
      end
      if ~isequal(local_data.iMoMo.PRCP,data.iMoMo.PRCP)
        diff = local_data.iMoMo.PRCP-data.iMoMo.PRCP;
        if any(any(abs(diff) > 0.001))
          doublefprintf(fid,'\n   WARNING: PRCP of %s nod equal %s within tolerance.\n',local_file,file);
        end
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
  
  elseif strcmp('WMO.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'raw',filesep,'WMO.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(local_data.WMO.ID,data.WMO.ID)
        doublefprintf(fid,'\n   WARNING: ID of %s not equal %s.\n',local_file,file);
      end
      if ~isequal(local_data.WMO.TEMP,data.WMO.TEMP)
        abs_diff = local_data.WMO.TEMP-data.WMO.TEMP;
        rel_diff1 = abs_diff ./ local_data.WMO.TEMP;
        rel_diff2 = abs_diff ./ data.WMO.TEMP;
        if any(any(abs(rel_diff1) > 0.5 & abs(rel_diff2) > 0.5))
          doublefprintf(fid,'\n   WARNING: TEMP of %s nod equal %s within tolerance.\n',local_file,file);
          keyboard
        end
      end
      if ~isequal(local_data.WMO.PRCP,data.WMO.PRCP)
        diff = local_data.WMO.PRCP-data.WMO.PRCP;
        if any(any(abs(diff) > 0.001))
          doublefprintf(fid,'\n   WARNING: PRCP of %s nod equal %s within tolerance.\n',local_file,file);
        end
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
  
  % Testing processRaw.
  elseif strcmp('SI_P_F.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'SI_P_F.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(local_data.a_P_F,data.a_P_F)
        diff = local_data.a_P_F - data.a_P_F;
        doublefprintf(fid,'\n   WARNING: a_P_F not equal. Testing tolerance.\n');
        if any(any(abs(diff) > 100))
          doublefprintf(fid,'\n   WARNING: a_P_F of %s not equal %s within tolerance.\n',local_file,file);
        end
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
    
  elseif strcmp('SI_P.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'SI_P.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(size(local_data.a_P),size(data.a_P))
        doublefprintf(fid,'\n   WARNING: size(a_P) of % not equal %s.\n',local_file,file);
      end
      for j = size(local_data.a_P,2):-2:1
        if ~isequal(local_data.a_P{j},data.a_P{j})
          diff = local_data.a_P{j} - data.a_P{j};
          doublefprintf(fid,'\n   WARNING: a_P{%d} not equal. Testing tolerance.\n',j);
          if any(any(abs(diff) > 100))
            doublefprintf(fid,'\n   WARNING: a_P{%d} of %s not equal %s within tolerance.\n',j,local_file,file);
          end
        end
      end
      
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
    
  elseif strcmp('SI_pET.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'SI_pET.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(size(local_data.a_pET),size(data.a_pET))
        doublefprintf(fid,'\n   WARNING: size(a_pET) of % not equal %s.\n',local_file,file);
      end
      for j = size(local_data.a_pET,2):-2:1
        if ~isequal(local_data.a_pET{j},data.a_pET{j})
          diff = local_data.a_pET{j} - data.a_pET{j};
          doublefprintf(fid,'\n   WARNING: a_pET{%d} not equal. Testing tolerance.\n',j);
          if any(any(abs(diff) > 100))
            doublefprintf(fid,'\n   WARNING: a_pET{%d} of %s not equal %s within tolerance.\n',j,local_file,file);
          end
        end
      end
      
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
    
  elseif strcmp('SI_pET_F.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'SI_pET_F.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(local_data.a_pET_F,data.a_pET_F)
        diff = local_data.a_pET_F - data.a_pET_F;
        doublefprintf(fid,'\n   WARNING: a_pET_F not equal. Testing tolerance.\n');
        if any(any(abs(diff) > 100))
          doublefprintf(fid,'\n   WARNING: a_pET_F of %s not equal %s within tolerance.\n',local_file,file);
        end
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
    
  elseif strcmp('SI_T.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'SI_T.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(size(local_data.a_T),size(data.a_T))
        doublefprintf(fid,'\n   WARNING: size(a_T) of % not equal %s.\n',local_file,file);
      end
      for j = size(local_data.a_T,2):-2:1
        if ~isequal(local_data.a_T{j},data.a_T{j})
          diff = local_data.a_T{j} - data.a_T{j};
        doublefprintf(fid,'\n   WARNING: a_T{%d} not equal. Testing tolerance.\n',j);
          if any(any(abs(diff) > 100))
            doublefprintf(fid,'\n   WARNING: a_T{%d} of %s not equal %s within tolerance.\n',j,local_file,file);
          end
        end
      end
      
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end

  elseif strcmp('sub_FP.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'sub_FP.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(size(local_data.FPsub),size(data.FPsub))
        doublefprintf(fid,'\n   WARNING: size(FPsub) of % not equal %s.\n',local_file,file);
      end
      for j = size(local_data.FPsub,2):-2:1
        if any(any(any(isnan(local_data.FPsub{j}.FP)))) || any(any(any(isnan(data.FPsub{j}.FP))))
          continue
        end
        if ~isequal(local_data.FPsub{j}.FP,data.FPsub{j}.FP)
          diff = local_data.FPsub{j}.FP - data.FPsub{j}.FP;
          doublefprintf(fid,'\n   WARNING: FPsub{%d}.FP not equal. Testing tolerance.\n',j);
          keyboard
          if any(any(abs(diff) > 100))
            doublefprintf(fid,'\n   WARNING: FPsub{%d}.FP of %s not equal %s within tolerance.\n',j,local_file,file);
          end
        end
        if ~isequal(local_data.FPsub{j}.timeF,data.FPsub{j}.timeF)
          doublefprintf(fid,'\n   WARNING FPsub{%d}.timeF of %s not equal %s.\n',j,local_file,file);
        end
      end
      
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
    
  elseif strcmp('sub_FpET.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'sub_FpET.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(size(local_data.FpETsub),size(data.FpETsub))
        doublefprintf(fid,'\n   WARNING: size(FpETsub) of % not equal %s.\n',local_file,file);
      end
      for j = size(local_data.FpETsub,2):-2:1
        if any(any(any(isnan(local_data.FpETsub{j}.FpET)))) || any(any(any(isnan(data.FpETsub{j}.FpET))))
          continue
        end
        if ~isequal(local_data.FpETsub{j}.FpET,data.FpETsub{j}.FpET)
          diff = local_data.FpETsub{j}.FpET - data.FpETsub{j}.FpET;
          doublefprintf(fid,'\n   WARNING: FpETsub{%d}.FpET not equal. Testing tolerance.\n',j);
          if any(any(abs(diff) > 100))
            doublefprintf(fid,'\n   WARNING: FpETsub{%d}.FpET of %s not equal %s within tolerance.\n',j,local_file,file);
          end
        end
        if ~isequal(local_data.FpETsub{j}.timeF,data.FpETsub{j}.timeF)
          doublefprintf(fid,'\n   WARNING FpETsub{%d}.timeF of %s not equal %s.\n',j,local_file,file);
        end
      end
      
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end

  elseif strcmp('sub_FT.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'sub_FT.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequal(size(local_data.FTsub),size(data.FTsub))
        doublefprintf(fid,'\n   WARNING: size(FTsub) of % not equal %s.\n',local_file,file);
      end
      for j = size(local_data.FTsub,2):-2:1
        if any(any(any(isnan(local_data.FTsub{j}.FT)))) || any(any(any(isnan(data.FTsub{j}.FT))))
          continue
        end
        if ~isequal(local_data.FTsub{j}.FT,data.FTsub{j}.FT)
          diff = local_data.FTsub{j}.FT - data.FTsub{j}.FT;
          doublefprintf(fid,'\n   WARNING: sub_FTsub{%d}.FT not equal. Testing tolerance.\n',j);
          if any(any(abs(diff) > 100))
            doublefprintf(fid,'\n   WARNING: FTsub{%d}.FT of %s not equal %s within tolerance.\n',j,local_file,file);
          end
        end
        if ~isequal(local_data.FTsub{j}.timeF,data.FTsub{j}.timeF)
          doublefprintf(fid,'\n   WARNING FTsub{%d}.timeF of %s not equal %s.\n',j,local_file,file);
        end
      end
      
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end

  elseif strcmp('sub_P.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'sub_P.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequaln(local_data.Psub,data.Psub)
        diff = local_data.Psub - data.Psub;
        doublefprintf(fid,'\n   WARNING: sub_Psub not equal. Testing tolerance.\n');
        keyboard
        if any(any(abs(diff) > 100))
          doublefprintf(fid,'\n   WARNING: Psub of %s not equal %s within tolerance.\n',local_file,file);
        end
      end
      if ~isequal(local_data.timeP,data.timeP)
        doublefprintf(fid,'\n   WARNING: timeP of %s not equal %s.\n',local_file,file);
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
 
  elseif strcmp('sub_pET.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'sub_pET.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequaln(local_data.pETsub,data.pETsub)
        diff = local_data.pETsub - data.pETsub;
        doublefprintf(fid,'\n   WARNING: sub_pETsub not equal. Testing tolerance.\n');
        if any(any(abs(diff) > 100))
          doublefprintf(fid,'\n   WARNING: pETsub of %s not equal %s within tolerance.\n',local_file,file);
        end
      end
      if ~isequal(local_data.timeT,data.timeT)
        doublefprintf(fid,'\n   WARNING: timeP of %s not equal %s.\n',local_file,file);
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end
    
  elseif strcmp('sub_T.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'sub_T.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequaln(local_data.sub_T,data.sub_T)
        diff = local_data.sub_T - data.sub_T;
        doublefprintf(fid,'\n   WARNING: sub_T not equal. Testing tolerance.\n');
        if any(any(abs(diff) > 100))
          doublefprintf(fid,'\n   WARNING: sub_T of %s not equal %s within tolerance.\n',local_file,file);
        end
      end
      if ~isequal(local_data.timeT,data.timeT)
        doublefprintf(fid,'\n   WARNING: timeT of %s not equal %s.\n',local_file,file);
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end

  elseif strcmp('sub_maxT.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'sub_maxT.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequaln(local_data.sub_maxT,data.sub_maxT)
        diff = local_data.sub_maxT - data.sub_maxT;
        doublefprintf(fid,'\n   WARNING: sub_maxT not equal. Testing tolerance.\n');
        if any(any(abs(diff) > 100))
          doublefprintf(fid,'\n   WARNING: sub_maxT of %s not equal %s within tolerance.\n',local_file,file);
        end
      end
      if ~isequal(local_data.timeT,data.timeT)
        doublefprintf(fid,'\n   WARNING: timeT of %s not equal %s.\n',local_file,file);
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end

  elseif strcmp('sub_minT.mat',file)
    local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
              'Themi',filesep,'data',filesep,'processed',filesep,'sub',filesep,...
              'sub_minT.mat'];
    if exist(local_file,'file') == 2 % Is a file.
      %doublefprintf(fid,'   File found in path.');
      % Load both files and compare them.
      local_data = load(local_file);
      data = load([test_dir,filesep,file]);
      if ~isequaln(local_data.sub_minT,data.sub_minT)
        diff = local_data.sub_minT - data.sub_minT;
        doublefprintf(fid,'\n   WARNING: sub_minT not equal. Testing tolerance.\n');
        if any(any(abs(diff) > 100))
          doublefprintf(fid,'\n   WARNING: sub_minT of %s not equal %s within tolerance.\n',local_file,file);
        end
      end
      if ~isequal(local_data.timeT,data.timeT)
        doublefprintf(fid,'\n   WARNING: timeT of %s not equal %s.\n',local_file,file);
      end
      doublefprintf(fid,'. . . done.\n');
    else
      doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
      return
    end

%   % Testing runModel.
%   elseif ~isempty(regexp(file,'^\d{6}\.mat$','once')) % Looks for a file with a name matching dddddd.mat.
%     local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
%               'Themi',filesep,'data',filesep,'processed',filesep,'db',filesep,...
%               file];
%     if exist(local_file,'file') == 2
%       local_data = load(local_file);
%       data = load([test_dir,filesep,file]);
%       if ~isequal(local_data.dateMet,data.dateMet)
%         doublefprintf(fid,'\n   WARNING: dateMet of %s not equal %s.\n',local_file,file);
%       end
%       if ~isequal(local_data.dateSim,data.dateSim)
%         doublefprintf(fid,'\n   WARNING: dateSim of %s not equal %s.\n',local_file,file);
%       end
%       if ~isequaln(local_data.ET,data.ET)
%         abs_diff = local_data.ET - data.ET;
%         rel_diff1 = abs_diff ./ local_data.ET;
%         rel_diff2 = abs_diff ./ data.ET;
%         doublefprintf(fid,'\n   WARNING: ET not equal. Testing tolerance.\n');
%         if any(any(any(abs(abs_diff) > 0.04 & abs(rel_diff1) > 0.05 & abs(rel_diff2) > 0.05)))
%           doublefprintf(fid,'\n   WARNING: ET of %s not equal %s within tolerance.\n',local_file,file);
% %           keyboard
%         end
%       end
%       if ~isequaln(local_data.G,data.G)
%         local_data.G(local_data.G < 0.0001) = 0;
%         data.G(data.G < 0.0001) = 0;
%         abs_diff = local_data.G - data.G;
%         rel_diff1 = abs_diff ./ local_data.G;
%         rel_diff2 = abs_diff ./ data.G;
%         doublefprintf(fid,'\n   WARNING: G not equal. Testing tolerance.\n');
%         if any(any(any(abs(abs_diff) > 0.04 & abs(rel_diff1) > 0.05 & abs(rel_diff2) > 0.05)))
%           doublefprintf(fid,'\n   WARNING: G of %s not equal %s within tolerance.\n',local_file,file);
% %           keyboard
%         end
%       end
%       if ~isequaln(local_data.T,data.T)
%         abs_diff = local_data.T - data.T;
%         rel_diff1 = abs_diff ./ local_data.T;
%         rel_diff2 = abs_diff ./ data.T;
%         doublefprintf(fid,'\n   WARNING: T not equal. Testing tolerance.\n');
%         if any(any(any(abs(abs_diff) > 0.04 & abs(rel_diff1) > 0.05 & abs(rel_diff2) > 0.05)))
%           doublefprintf(fid,'\n   WARNING: T of %s not equal %s within tolerance.\n',local_file,file);
% %           keyboard
%         end
%       end
%       if ~isequaln(local_data.P,data.P)
%         abs_diff = local_data.P - data.P;
%         rel_diff1 = abs_diff ./ local_data.P;
%         rel_diff2 = abs_diff ./ data.P;
%         doublefprintf(fid,'\n   WARNING: P not equal. Testing tolerance.\n');
%         if any(any(any(abs(abs_diff) > 0.04 & abs(rel_diff1) > 0.05 & abs(rel_diff2) > 0.05)))
%           doublefprintf(fid,'\n   WARNING: P of %s not equal %s within tolerance.\n',local_file,file);
% %           keyboard
%         end
%       end
%       if ~isequaln(local_data.Q,data.Q)
%         local_data.Q(local_data.Q < 0.0001) = 0;
%         data.Q(data.Q < 0.0001) = 0;
%         abs_diff = local_data.Q - data.Q;
%         rel_diff1 = abs_diff ./ local_data.Q;
%         rel_diff2 = abs_diff ./ data.Q;
%         doublefprintf(fid,'\n   WARNING: Q not equal. Testing tolerance.\n');
%         if any(any(any( abs(abs_diff) > 0.04 & abs(rel_diff1) > 0.05 & abs(rel_diff2) > 0.05 )))
%           doublefprintf(fid,'\n   WARNING: Q of %s not equal %s within tolerance.\n',local_file,file);
% %           keyboard
%         end
%       end
%       if ~isequaln(local_data.S,data.S)
%         abs_diff = local_data.S - data.S;
%         rel_diff1 = abs_diff ./ local_data.S;
%         rel_diff2 = abs_diff ./ data.S;
%         doublefprintf(fid,'\n   WARNING: S not equal. Testing tolerance.\n');
%         if any(any(any(abs(abs_diff) > 0.04 & abs(rel_diff1) > 0.05 & abs(rel_diff2) > 0.05)))
%           doublefprintf(fid,'\n   WARNING: S of %s not equal %s within tolerance.\n',local_file,file);
% %           keyboard
%         end
%       end
%       doublefprintf(fid,'. . . done.\n');
%     else
%       doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
%       return
%     end
%     
%   elseif ~isempty(regexp(file,'^\d{6}_E\.mat$','once')) % Looks for a file with a name matching dddddd.mat.
%     local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
%               'Themi',filesep,'results',filesep,'nowcast',filesep,file];
%     if exist(local_file,'file') == 2
%       local_data = load(local_file);
%       data = load([test_dir,filesep,file]);
%       if ~isequal(local_data.E,data.E)
%         doublefprintf(fid,'\n   WARNING: E not equal. Testing tolerance.\n');
%         local_data.E(local_data.E < 0.001) = 0;
%         data.E(data.E < 0.001) = 0;
%         abs_diff = local_data.E - data.E;
%         rel_diff = abs_diff ./ local_data.E;
%         rel_diff2 = abs_diff ./ data.E;
%         rel_diff(isnan(rel_diff)) = 0;
%         rel_diff2(isnan(rel_diff2)) = 0;
%         if any(any(abs(abs_diff) > 0.04 & abs(rel_diff) > 0.05 & abs(rel_diff2) > 0.05))
%           doublefprintf(fid,'\n   WARNING: E of %s not equal %s within tolerance.\n',local_file,file);
% %           keyboard
%         end
%       end
%       doublefprintf(fid,'. . . done.\n');
%     else
%       doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
%       return
%     end
%       
%   elseif ~isempty(regexp(file,'^\d{6}_FET\.mat$','once')) % Looks for a file with a name matching dddddd.mat.
%     local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
%               'Themi',filesep,'results',filesep,'forecast',filesep,'ET',...
%               filesep,file];
%     if exist(local_file,'file') == 2
%       local_data = load(local_file);
%       data = load([test_dir,filesep,file]);
%       if ~isequal(local_data.FET.timeF,data.FET.timeF)
%         doublefprintf(fid,'\n   WARNING FET.timeF of %s not equal %s.\n', local_file,file);
%       end
%       if ~isequal(size(local_data.FET.catch,2),size(data.FET.catch,2))
%         doublefprintf(fid,'   ERROR: size(FET.catch).\n');
%         return
%       else
%         for j = size(local_data.FET.catch,2):-5:1
%           if ~isequal(size(local_data.FET.catch{j}.Q,2),size(data.FET.catch{j}.Q,2))
%             doublefprintf(fid,'   ERROR: size(FET.catch{%d}.Q).\n',j);
%             return
%           else
%             for k = size(local_data.FET.catch{j}.Q,2):-5:1
%               if ~isequal(local_data.FET.catch{j}.Q{k},data.FET.catch{j}.Q{k})
%                 doublefprintf(fid,'\n   WARNING: FET.catch{%d}.Q{%d} not equal. Testing tolerance.\n',j,k);
%                 abs_diff = abs(local_data.FET.catch{j}.Q{k} - data.FET.catch{j}.Q{k});
%                 rel_diff1 = abs_diff ./ abs(local_data.FET.catch{j}.Q{k});
%                 rel_diff2 = abs_diff ./ abs(data.FET.catch{j}.Q{k});
%                 if any(any(abs_diff > 0.04 & rel_diff1 > 0.05 & rel_diff2 > 0.05))
%                   doublefprintf(fid,'\n   WARNING: FET.catch{%d}.Q{%d} of %s not equal %s within tolerance.\n',j,k,local_file,file);
% %                   keyboard
%                 end
%               end
%             end
%           end
%         end
%       end
%       doublefprintf(fid,'. . . done.\n');
%     else
%       doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
%       return
%     end
%     
%   elseif ~isempty(regexp(file,'^\d{6}_FG\.mat$','once')) % Looks for a file with a name matching dddddd.mat.
%     local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
%               'Themi',filesep,'results',filesep,'forecast',filesep,'G',...
%               filesep,file];
%     if exist(local_file,'file') == 2
%       local_data = load(local_file);
%       data = load([test_dir,filesep,file]);
%       if ~isequal(local_data.FG.timeF,data.FG.timeF)
%         doublefprintf(fid,'\n   WARNING FG.timeF of %s not equal %s.\n', local_file,file);
%       end
%       if ~isequal(size(local_data.FG.catch,2),size(data.FG.catch,2))
%         doublefprintf(fid,'   ERROR: size(FG.catch).\n');
%         return
%       else
%         for j = size(local_data.FG.catch,2):-5:1
%           if ~isequal(size(local_data.FG.catch{j}.Q,2),size(data.FG.catch{j}.Q,2))
%             doublefprintf(fid,'   ERROR: size(FG.catch{%d}.Q).\n',j);
%             return
%           else
%             for k = size(local_data.FG.catch{j}.Q,2):-5:1
%               if ~isequal(local_data.FG.catch{j}.Q{k},data.FG.catch{j}.Q{k})
%                 diff = local_data.FG.catch{j}.Q{k} - data.FG.catch{j}.Q{k};
%                 if any(any(abs(diff) > 0.6))
%                   doublefprintf(fid,'\n   WARNING: FG.catch{%d}.Q{%d} of %s not equal %s within tolerance.\n',j,k,local_file,file);
% %                   keyboard
%                 end
%               end
%             end
%           end
%         end
%       end
%       doublefprintf(fid,'. . . done.\n');
%     else
%       doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
%       return
%     end
% 
%   elseif ~isempty(regexp(file,'^\d{6}_FQ\.mat$','once')) % Looks for a file with a name matching dddddd.mat.
%     local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
%               'Themi',filesep,'results',filesep,'forecast',filesep,'Q',...
%               filesep,file];
%     if exist(local_file,'file') == 2
%       local_data = load(local_file);
%       data = load([test_dir,filesep,file]);
%       if ~isequal(local_data.FQ.timeF,data.FQ.timeF)
%         doublefprintf(fid,'\n   WARNING FQ.timeF of %s not equal %s.\n', local_file,file);
%       end
%       if ~isequal(size(local_data.FQ.catch,2),size(data.FQ.catch,2))
%         doublefprintf(fid,'   ERROR: size(FQ.catch).\n');
%         return
%       else
%         for j = size(local_data.FQ.catch,2):-5:1
%           if ~isequal(size(local_data.FQ.catch{j}.Q,2),size(data.FQ.catch{j}.Q,2))
%             doublefprintf(fid,'   ERROR: size(FQ.catch{%d}.Q).\n',j);
%             return
%           else
%             for k = size(local_data.FQ.catch{j}.Q,2):-5:1
%               if ~isequal(local_data.FQ.catch{j}.Q{k},data.FQ.catch{j}.Q{k})
%                 diff = local_data.FQ.catch{j}.Q{k} - data.FQ.catch{j}.Q{k};
%                 if any(any(abs(diff) > 0.6))
%                   doublefprintf(fid,'\n   WARNING: FQ.catch{%d}.Q{%d} of %s not equal %s within tolerance.\n',j,k,local_file,file);
% %                   keyboard
%                 end
%               end
%             end
%           end
%         end
%       end
%       doublefprintf(fid,'. . . done.\n');
%     else
%       doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
%       return
%     end
%     
%   elseif ~isempty(regexp(file,'^\d{6}_FS\.mat$','once')) % Looks for a file with a name matching dddddd.mat.
%     local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
%               'Themi',filesep,'results',filesep,'forecast',filesep,'S',...
%               filesep,file];
%     if exist(local_file,'file') == 2
%       local_data = load(local_file);
%       data = load([test_dir,filesep,file]);
%       if ~isequal(local_data.FS.timeF,data.FS.timeF)
%         doublefprintf(fid,'\n   WARNING FS.timeF of %s not equal %s.\n', local_file,file);
%       end
%       if ~isequal(size(local_data.FS.catch,2),size(data.FS.catch,2))
%         doublefprintf(fid,'   ERROR: size(FS.catch).\n');
%         return
%       else
%         for j = size(local_data.FS.catch,2):-5:1
%           if ~isequal(size(local_data.FS.catch{j}.Q,2),size(data.FS.catch{j}.Q,2))
%             doublefprintf(fid,'   ERROR: size(FS.catch{%d}.Q).\n',j);
%             return
%           else
%             for k = size(local_data.FS.catch{j}.Q,2):-5:1
%               if ~isequal(local_data.FS.catch{j}.Q{k},data.FS.catch{j}.Q{k})
%                 diff = local_data.FS.catch{j}.Q{k} - data.FS.catch{j}.Q{k};
%                 if any(any(abs(diff) > 0.6))
%                   doublefprintf(fid,'\n   WARNING: FS.catch{%d}.Q{%d} of %s not equal %s within tolerance.\n',j,k,local_file,file);
%                   keyboard
%                 end
%               end
%             end
%           end
%         end
%       end
%       doublefprintf(fid,'. . . done.\n');
%     else
%       doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
%       return
%     end
%   
%   elseif strcmp(file,'E.mat') % check if initial states are the same.
%     local_file = [current_dir,filesep,'..',filesep,'..',filesep,'app',filesep,...
%               'Themi',filesep,'resources',filesep,'restart',filesep,'E.mat'];
%     if exist(local_file,'file') == 2
%       local_data = load(local_file);
%       data = load([test_dir,filesep,file]);
%       if ~isequaln(local_data.E,data.E)
%         abs_diff = abs( local_data.E - data.E );
%         if any(any(abs_diff > 0.001))
%           doublefprintf(fid,'\n   WARNING: Initial state files are not the same! Continuing.\n');
% %           keyboard
% %           return
%         end
%       end
%       doublefprintf(fid,'. . . done.\n');
%     else
%       doublefprintf(fid,'   ERROR: File %s not found in local path %s.\n',file,local_file);
%       return
%     end
    
  else
    doublefprintf(fid,'ERROR file %s not found.\n',file);
  end
  
end

fclose(fid);
end % compare_data



%% Local functions

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