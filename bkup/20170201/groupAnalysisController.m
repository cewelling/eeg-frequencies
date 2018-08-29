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

%% Group SNR Analysis------------------------------------------------------

% Has the participant undergone SNR Analysis?

if snrAnal == 1
    compareSNRs(groupParNums) 
end
SNRs = {};
for iGroup = 1:length(groupParNums)
    group = groupParNums{iGroup};
    SNRs{iGroup} = [];
    for parNum = group
        for currentRun = runIndices
            % Ensure that SNR has been calculated for this participant/electrode set/run
%             if ~exist(['indices/SNR/stratus' num2str(parNum) '_run' num2str(currentRun) '_' electrodeSet.name '.mat'])
%                 error(['SNR has not been calculated for stratus' num2str(parNum) ' run ' num2str(currentRun) ' with ' electrodeSet.name])
%             end
            
            % load and store data
            if exist(['indices/SNR/stratus' num2str(parNum) '_run' num2str(currentRun) '_' electrodeSet.name '.mat'])
                load (['indices/SNR/stratus' num2str(parNum) '_run' num2str(currentRun) '_' electrodeSet.name '.mat'])
            else
                continue;
            end
            
            for i = 1:length(meanSNRs)
                
                % Choose SNRs calculated with the desired frequency and noise window
                if (meanSNRs(i).freq == FOI && meanSNRs(i).noiseWindowHalf == noiseWindowHalf)
                    if (meanSNRs(i).value > minSNR)
                        SNRs{iGroup} = [SNRs{iGroup} meanSNRs(i).value];
                    end
                end
            end
        end
    end
end

% Calculate average SNR and print results
formatSpec = '%s ''s mean SNR Index is %4.3f\n';
formatSpec2 = '%d runs included\n\n';
for iGroup = 1:length(stratusGroups)
    meanSNR = mean(SNRs{iGroup});   
    fprintf(formatSpec, stratusGroupIDs{iGroup}, meanSNR)  
    fprintf(formatSpec2, length(SNRs{iGroup}))
end

H = ttest2(SNRs{1}, SNRs{2});
if H ~= 0
    disp('SNRs are significantly different')
else
    disp('SNRs are not significantly different')
end
 
%% Group plotting

makeGroupPlot();

%% Statistics

% indexAnal

%%
