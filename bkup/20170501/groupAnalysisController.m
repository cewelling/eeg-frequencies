% Group analysis of binocular rivalry dynamics
% aggregates and analyzes results produced by analysisController.m

clearvars

%% Switches----------------------------------------------------------------

% load parameters
groupAnalysisParams

%% Prep--------------------------------------------------------------------

% Weed out participants that don't exist (allows us to conveniently input
% full ranges for stratusGroups and cumulusGroups)
parsedGroup = [];
% for iGroup = 1:length(groupParNums)
%     parNums = groupParNums{iGroup};
%     parsedGroup = [];
%     for parNum = parNums
%         if exist(['pilotlegends/' groupCodes{iGroup}' num2str(parNum) '.txt'],'file')
%             parsedGroup = [parsedGroup parNum];
%         end
%     end
%     groupParNums{iGroup} = parsedGroup;
% end
    
%% Check whether a participant has undergone all of the necessary analysis
% steps in analysisController.m 
% Print steps that still need to be done...actually this can be done
% throughout the script

%% Gate Runs---------------------------------------------------------------

runsIncluded = gateRuns();
% return a list of the participants and runs that were included in the
% analysis


 
%% Group plotting

makeGroupPlot();

%% Statistics

% indexAnal

%%
