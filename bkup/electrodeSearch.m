clear all; close all; clc;
% setFigProps

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveIt = 'no';

% pars = [105;106;109;114;115];
pars = [];

pars(1).ID = 105;
pars(1).runs = [7 8];
pars(2).ID = 106;
pars(2).runs = [7 8];
pars(3).ID = 109;
pars(3).runs = [7 8 9];
pars(4).ID = 114;
pars(4).runs = [7 8 9];
pars(5).ID = 118;
pars(5).runs = [7 8 9];
% pars(5).ID = 115;
% pars(5).runs = [5];

envelopes = [];
envelopes(1).ID = 105;
envelopes(2).ID = 106;
envelopes(3).ID = 109;
envelopes(4).ID = 114;
envelopes(5).ID = 118;

electsOfInt = [5:16,21:24];

% add electrode fields to par struct

for iPar = 1:size(pars,2)
    for iElect = 1:length(electsOfInt)
        electNum = ['elect' num2str(electsOfInt(iElect))];
        pars(iPar).(electNum) = [0 0];
        envelopes(iPar).(electNum) = [];
    end
end

% runsOfInt = [7 8 9]; %[7 8 11 12]; %[1 4 5 6];
runType = 'rival'; % 'rival' or 'sim'
stimType = 'marz'; % 'darts' or 'marz'

epochHalf = 2; %+/- from the time of the button press %1.5
winLengthRLS = 750; %time in ms %1000
normType = 'z'; %'z' or 'norm' or 'none'
% analElect = 29; % use map to select electrode of interest. Default, 'Oz' = 29'

for parIdx = 1:size(pars,2)
    par = pars(parIdx).ID;
    runsOfInt = pars(parIdx).runs;

for analElectIdx = 1:length(electsOfInt) %change cab input to reflect this
    analElect = electsOfInt(analElectIdx);

loopCount = 0;
    
for analRunNum = runsOfInt %[7 8 11 12] % [1 4 5 6]

%%%%% Set up: Things to change %%%%%%%%%%%%%%%
% analRunNum = 6;
parName = ['stratus' num2str(par)];
fileNames = dir(['../../runScripts/Results/' parName '*']);
legend = importdata(['pilotlegends/' parName '.txt']); %pilot-' date '-caroline.txt
% stimulation_frequencies = [14.16 17]; %[17 21.25] Note: this is redefined later....
intermodulation_frequencies = 11.3;

%%%%% Set up: Always done %%%%%%%%%%%%%%%%%%%%%
fileNames = {fileNames.name}; filenames{1}(1).name = char(fileNames(analRunNum));
stimulation_frequencies = legend.data(analRunNum,1:2);
numTrials = legend.data(analRunNum,3);
trialDur = legend.data(analRunNum,4); 
defaultAnalParams %Load default anal params

stimulation_frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[maxResponseElecs, maxIndices] = runFFT(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend)

omniData = runKeyAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,minEpoch);

if strcmp(runType,'rival')
    if loopCount == 0
        [slowDurs, fastDurs, mixedDurs] = pressDurations(omniData,parName,analRunNum);
    else
        try
            [newslowDurs, newfastDurs, newmixedDurs] = pressDurations(omniData,parName,analRunNum);
        catch
            [slowDurs, fastDurs] = pressDurations(omniData,parName,analRunNum);
        end
    end
elseif strcmp(runType,'sim')
    if loopCount == 0
        [slowDurs, fastDurs] = pressDurations(omniData,parName,analRunNum);
    else
        [newslowDurs, newfastDurs] = pressDurations(omniData,parName,analRunNum);
    end
end

[multiFreqSeries, trialTime] = compareRLSfast(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend, omniData,winLengthRLS,analElect);


if loopCount == 0
    [transitions] = compareEpochs(omniData, multiFreqSeries, trialTime,analRunNum,epochHalf,normType);
else
    [newTransitions] = compareEpochs(omniData, multiFreqSeries, trialTime,analRunNum,epochHalf,normType);
end
    
%amend the first run data with new
    
if loopCount ~= 0
    for i = 1:6
        try
            transitions(i).f1 = [transitions(i).f1; newTransitions(i).f1];
            transitions(i).f2 = [transitions(i).f2; newTransitions(i).f2];
        catch
            transitions(i).f1 = transitions(i).f1;
            transitions(i).f2 = transitions(i).f2;
        end
            
    end
end

if strcmp(runType,'rival')
    if loopCount ~= 0
        slowDurs = [slowDurs; newslowDurs];
        fastDurs = [fastDurs; newfastDurs];
        mixedDurs = [mixedDurs; newmixedDurs];
    end
elseif strcmp(runType,'sim')
    if loopCount ~= 0
        slowDurs = [slowDurs; newslowDurs];
        fastDurs = [fastDurs; newfastDurs];
    end
end


loopCount = loopCount+1;

close all

end %finished aggrigating runs for that electrode

% determine if the transition signature is present

for iCase = 1:size(transitions,2)
    
        meanF1trace = mean(transitions(iCase).f1,1);
        errorF1trace = ste(transitions(iCase).f1);

        meanF2trace = mean(transitions(iCase).f2,1);
        errorF2trace = ste(transitions(iCase).f2);
        
        transitionTP = length(meanF1trace)/2;
        
        if iCase == 1
            if meanF2trace(transitionTP) > meanF1trace(transitionTP)
                % low to high checks
                pars(parIdx).(['elect' num2str(analElect)])(1) = 1;
            else
                % low to high does not check
                pars(parIdx).(['elect' num2str(analElect)])(1) = 0;
            end
        elseif iCase == 2
            if meanF1trace(transitionTP) > meanF2trace(transitionTP)
                % high to low checks
                pars(parIdx).(['elect' num2str(analElect)])(2) = 1;
            else
                % high to low does not check
                pars(parIdx).(['elect' num2str(analElect)])(2) = 0;
            end
        end
end

% strore the frequency envelopes for plotting later
if sum(pars(parIdx).(['elect' num2str(analElect)])) == 2
    envelopes(parIdx).(['elect' num2str(analElect)]) = transitions;
end
                
clearvars -except saveIt pars envelopes electsOfInt runType stimType epochHalf winLengthRLS normType runsOfInt par analElect parIdx             

end

end

%%% Plotting %%%

% close all

%RLS transition data

figure(analElect)
for iCase = 1:size(transitions,2)
    
        meanF1trace = mean(transitions(iCase).f1,1);
        errorF1trace = ste(transitions(iCase).f1);

        meanF2trace = mean(transitions(iCase).f2,1);
        errorF2trace = ste(transitions(iCase).f2);
    
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
            mseb([(-epochHalf+(1/512)):(1/512):epochHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
            vline(0,'k','transition')
            xlabel('Time from button press (s)');
            ylabel('Normalized amplitude');
        end
end
suptitle(['Stratus' num2str(pars) ' Transitions at Electrode ' num2str(analElect) ': ' runType '-' stimType])

saveName = ['Stratus' num2str(pars) '_Transitions-' runType '-' stimType];
saveas(gcf,['plots/transitions/' saveName],'pdf')


%Button press data
slowDur = mean(slowDurs);
slowError = std(slowDurs);

fastDur = mean(fastDurs);
fastError = std(fastDurs);

if strcmp(runType,'rival')
    mixedDur = mean(mixedDurs);
    mixedError = std(mixedDurs);

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
