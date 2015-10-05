function [] = testRunModel()
%% function [] = testRunModel()
%
% Calls different tests for runModel_Themi.m.
%


[bol,message] = testReadingSubAndForecastDataAndWritingDB();
if bol ~= 1
  fprintf('ERROR in testReadingSubAndForecastDataAndWritingDB().\n');
  fprintf('Message: %s\n',message);
end;

[bol, message] = testAssimilation();
if bol ~= 1
  fprintf('ERROR in testAssimilation().\n');
  fprintf('Message: %s.\n',message);
end

end