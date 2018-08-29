%function transitionAnal(pars, runsOfInt, runType, stimType, saveIt)

clearvars -except pars runsOfInt runType stimType saveIt; close all; clc;
% setFigProps

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveIt = 'no'; % 'yes' or 'no', toggle to output the .mat file for group anal

pars = [136]; % participant ID num - select which participant to run

runsOfInt = [4 5 6];  % select the runs (FROM LEGEND) to look at, should be from a single stim/trial type %rival/darts:[1 2 3], sim/darts:[4 5 6], rival/marz:[7 8 9]
runType = 'sim'; % 'rival' or 'sim', change this to analyze correclty and safe with appropriate title
stimType = 'darts'; % 'darts' or 'marz'

epochHalf = 5; %+/- from the time of the button press %1.5
winLengthRLS = 1000; %time in ms %1000
normType = 'z'; %'z' or 'norm' or 'none' %default = z-scoring ('z')
analElect = 29; % use map to select electrode of interest. Default, 'Oz' = 29'

plotWindow = 3; %amount of time to show in plots 
clearPlots = 'no'; % close out the plots generated - select 'no' to preserve the output of trial timecourse plots

loopCount = 0;
for par = pars
    
for analRunNum = runsOfInt


parName = ['stratus' num2str(par)];
fileNames = dir(['../../runScripts/Results/' parName '*']);
legend = importdata(['pilotlegends/' parName '.txt']); %pilot-' date '-caroline.txt
% stimulation_frequencies = [5.7 8.5]; %[17 21.25] Note: this is redefined later....
% intermodulation_frequencies = 11.3; This is defined in line 45, no?

fileNames = {fileNames.name}; filenames{1}(1).name = char(fileNames(analRunNum));
stimulation_frequencies = legend.data(analRunNum,1:2);
numTrials = legend.data(analRunNum,3);
trialDur = legend.data(analRunNum,4); 
defaultAnalParams %Load default anal params

intermodulation_frequencies = [(2*stimulation_frequencies(1) - stimulation_frequencies(2)), (2*stimulation_frequencies(2) - stimulation_frequencies(1))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run FFT - look at presence of stim frequencies in the spectrum
runFFT(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend)

% Run key anal - get key press data for transition analysis
[omniData,omniPresent] = runKeyAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,minEpoch);

% Get epoch durations
if loopCount == 0
    [slowDurs, fastDurs, mixedDurs] = pressDurations(omniData, parName, analRunNum);
else
    [newslowDurs, newfastDurs, newmixedDurs] = pressDurations(omniData, parName, analRunNum);
end

% I don't think we should have to define that the sim runs can't have mixed
% percepts...that should come out of the data 
% if strcmp(runType,'rival') % get epoch durations seperately for rival and sim runs - no mixed epochs for for simulation
%     if loopCount == 0
%         try
%             [slowDurs, fastDurs, mixedDurs] = pressDurations(omniData,parName,analRunNum);
%         catch
%             [slowDurs, fastDurs] = pressDurations(omniData,parName,analRunNum);
%             mixedDurs = [];
%         end
%     else
%         try
%             [newslowDurs, newfastDurs, newmixedDurs] = pressDurations(omniData,parName,analRunNum);
%         catch   % handle subjects with no mixed epochs
%             [newslowDurs, newfastDurs] = pressDurations(omniData,parName,analRunNum);
%             newmixedDurs = [];
%         end
%     end
% elseif strcmp(runType,'sim')
%     if loopCount == 0
%         [slowDurs, fastDurs] = pressDurations(omniData,parName,analRunNum);
%     else
%         [newslowDurs, newfastDurs] = pressDurations(omniData,parName,analRunNum);
%     end
% end

[multiFreqSeries, trialTime] = compareRLSfast(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend, omniData,winLengthRLS,analElect);

% compareEpochs is accomplishing the transition analysis and plotting
if loopCount == 0
    [transitions] = compareEpochs(omniData, multiFreqSeries, trialTime,analRunNum,epochHalf,normType);
else
    [newTransitions] = compareEpochs(omniData, multiFreqSeries, trialTime,analRunNum,epochHalf,normType);
end
    
%amend the first run data with new
    
if loopCount ~= 0
    for i = 1:6 % number of transition types
        try
            transitions(i).f1 = [transitions(i).f1; newTransitions(i).f1];
            transitions(i).f2 = [transitions(i).f2; newTransitions(i).f2];
        catch
            transitions(i).f1 = transitions(i).f1; % if a type of transition doesn't yet have any data
            transitions(i).f2 = transitions(i).f2;
        end
            
    end
end


if loopCount ~= 0
    slowDurs = [slowDurs; newslowDurs];
    fastDurs = [fastDurs; newfastDurs];
    mixedDurs = [mixedDurs; newmixedDurs];
end

% if strcmp(runType,'rival')
%     if loopCount ~= 0
%         slowDurs = [slowDurs; newslowDurs];
%         fastDurs = [fastDurs; newfastDurs];
%         mixedDurs = [mixedDurs; newmixedDurs];
%     end
% elseif strcmp(runType,'sim')
%     if loopCount ~= 0
%         slowDurs = [slowDurs; newslowDurs];
%         fastDurs = [fastDurs; newfastDurs];
%     end
% end


loopCount = loopCount+1;

if strcmp(clearPlots,'yes');
    close all % clear out all the plots
end

end
end

%%% Plotting %%%

% close all

%RLS transition data

%settings
lineClr = ['r';'g'];
lineProps = [];
lineProps.col = cellstr(lineClr);

%loop for transition type

multiplot = figure;
for iCase = 1:size(transitions,2)
    
        meanF1trace = nanmean(transitions(iCase).f1,1);
        errorF1trace = ste(transitions(iCase).f1);

        meanF2trace = nanmean(transitions(iCase).f2,1);
        errorF2trace = ste(transitions(iCase).f2);
    
        figure(multiplot)
        if length(errorF1trace) > 1
            switch runType
                case 'rival'
                    if iCase < 2 || iCase > 5
                        subplot(2,3,iCase)
                    elseif iCase == 2 || iCase == 3
                        subplot(2,3,(iCase+2))
                    elseif iCase == 4 || iCase == 5
                        subplot(2,3,(iCase-2))
                    end
                case 'sim'
                    if iCase == 1
                        subplot(2,3,2)
                    elseif iCase == 2
                        subplot(2,3,5)
                    end
            end
            title(transitions(iCase).type)
            hold on
            mseb([(-epochHalf+(1/512)):(1/512):epochHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],lineProps,1);
            axis([-plotWindow plotWindow -1.5 1.5])
            vline(0,'k')
            xtick = get(gca,'XTick');
            set(gca,'FontSize',10)
            
            %Create single plot
            figure(iCase+1)
            title(transitions(iCase).type)
            hold on
            mseb([(-epochHalf+(1/512)):(1/512):epochHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],lineProps,1);
            axis([-plotWindow plotWindow -1.5 1.5])
            vline(0,'k') %(0,'k','transition')
            xtick = get(gca,'XTick');
            set(gca,'FontSize',10)
            xlabel('Time from Button Press (s)');
            ylabel('Frequency Band Power');
            
            %save single plot
            targDirectory = ['plots/transitions/Individual/' char(num2str(pars))];
            if ~exist(targDirectory, 'dir')
                mkdir('plots/transitions/Individual/', char(num2str(pars)));
            end
            saveName = char(strcat(num2str(pars), '/Stratus', num2str(pars), '-', transitions(iCase).type, '-', runType, '-', stimType));
            saveas(gcf,['plots/transitions/Individual/' saveName],'pdf')
        end
end

figure(multiplot)
suptitle(['Stratus' num2str(pars) ' Transitions at Electrode ' num2str(analElect) ': ' runType '-' stimType])
[ax,h1] = suplabel('Time from Button Press (s)','x',[.1 .1 .84 .84]); %[.08 .08 .84 .84]
[ay,h2] = suplabel('Frequency Band Power','y',[.15 .15 .6 .6]);
set(h1,'fontsize',15)
set(h2,'fontsize',15)

if strcmp('yes',saveIt)
    saveName = char(strcat(num2str(pars), '/Stratus', num2str(pars), '_Transitions-', runType, '-', stimType));
    saveas(gcf,['plots/transitions/Individual/' saveName],'pdf')
end

%Button press data
slowDur = nanmean(slowDurs);
slowError = nanstd(slowDurs);

fastDur = nanmean(fastDurs);
fastError = nanstd(fastDurs);

if strcmp(runType,'rival')
    mixedDur = nanmean(mixedDurs);
    mixedError = nanstd(mixedDurs);

    figure
    hold on
    bar(1:3,[slowDur fastDur mixedDur])
    Labels = {'5.67', '8.5 Hz', 'Mixed'};
    set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
    errorbar(1:3,[slowDur fastDur mixedDur],[slowError fastError mixedError],'.')
    title(['Average Percept Durations, ' parName])
    ylabel('seconds')
elseif strcmp(runType,'sim')
    figure
    hold on
    bar(1:2,[slowDur fastDur])
    Labels = {'5.67', '8.5 Hz'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    errorbar(1:2,[slowDur fastDur],[slowError fastError],'.')
    title(['Average Percept Durations, ' parName])
    ylabel('seconds')
end

if strcmp('yes',saveIt)
    saveName = char(strcat(num2str(pars), '/Stratus', num2str(pars), '_Durations-', runType, '-', stimType));
    saveas(gcf,['plots/transitions/Individual/' saveName],'pdf')
end

durations(1).type = 'Slow';
durations(1).duration = slowDur;
durations(2).type = 'Fast';
durations(2).duration = fastDur;
durations(3).type = 'Mixed';
if exist('mixedDur','var')
    durations(3).duration = mixedDur;
end

if strcmp('yes',saveIt)
    dataName = ['./epochResults/' parName '_' stimType '_' runType '_Transitions'];
    save(dataName,'transitions','durations');
end

% end