function [] = testRunModelOpenDA()
%% function [] = testRunModel()
%
% Calls different tests for runModelOpenDA_Themi.m.
%


[bol, message] = testAssimilationOpenDA();
if bol ~= 1
  fprintf('ERROR in testAssimilationOpenDA().\n');
  fprintf('Message: %s.\n',message);
end

end