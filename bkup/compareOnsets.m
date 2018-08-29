clear all; close all; clc;
% setFigProps

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runsOfInt = [7 8 11 12]; %[7 8 11 12]; %[1 4 5 6];
par = [100];


loopCount = 0;
for analRunNum = runsOfInt %[7 8 11 12] % [1 4 5 6]

%%%%% Set up: Things to change %%%%%%%%%%%%%%%
% analRunNum = 6;
parName = ['stratus' num2str(par)];
% date = '4_15';
% fileNames = dir(['../../runScripts/Results/' parName '*' date '*']);
fileNames = dir(['../../runScripts/Results/' parName '*']);
legend = importdata(['pilotlegends/' parName '.txt']); %pilot-' date '-caroline.txt
stimulation_frequencies = [14.16 17]; %[17 21.25] Note: this is redefined later....
intermodulation_frequencies = 11.3;

%%%%% Set up: Always done %%%%%%%%%%%%%%%%%%%%%
fileNames = {fileNames.name}; filenames{1}(1).name = char(fileNames(analRunNum));
stimulation_frequencies = legend.data(analRunNum,1:2);
numTrials = legend.data(analRunNum,3);
trialDur = legend.data(analRunNum,4); 
defaultAnalParams %Load default anal params

% stimulation_frequencies = [7.0833 8.5000];

stimulation_frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omniData = runKeyAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,minEpoch);

% sort the trials depending on which dominant state they enter first

if loopCount == 0
    onsets(1).type = 'low';
    onsets(2).type = 'high';
    lowStarts = [];
    highStarts = [];
end

thisRun = num2str(loopCount+1);
runName = ['run' thisRun];
onsets(1).(runName) = [];
onsets(2).(runName) = [];

for iTrial = 1:max(omniData(:,1))
    thesePercepts = omniData(find(omniData(:,1) == iTrial),:);
    perceptIndex = find(thesePercepts(:,7),1);
    perceptType = thesePercepts(perceptIndex,7);
    perceptStart = thesePercepts(perceptIndex,2);
    if perceptType < 0
        onsets(1).(runName) = [onsets(1).(runName), iTrial];
        lowStarts = [lowStarts, perceptStart];
    elseif perceptType > 0
        onsets(2).(runName) = [onsets(2).(runName), iTrial];
        highStarts = [highStarts, perceptStart];
    end
                 
end

multiFreqSeries = [];
trialTime = [];

[multiFreqSeries, trialTime] = compareRLSfast(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend, omniData);


%Normalize the frequency band time series
normalFreqSeries = [];

for iTrial = 1:size(multiFreqSeries,2)
    
    for iFreq = 1:size(multiFreqSeries{1,iTrial},2);
        
        RLSamp = cell2mat(multiFreqSeries{1,iTrial}(1,iFreq));
        
        normRLSamp = 2*(RLSamp - min(RLSamp))/(max(RLSamp) - min(RLSamp)) - 1;
        
        normalFreqSeries{1,iTrial}{1,iFreq} = normRLSamp;
        
    end
    
end

trialStartTime = 5; %how much time at the start of the trial to analyze

if loopCount == 0
    F1resolveLow = [];
    F2resolveLow = [];

    F1resolveHigh = [];
    F2resolveHigh = [];
end

for iTrial = 1:size(normalFreqSeries,2)
    
    RLSf1 = cell2mat(normalFreqSeries{1,iTrial}(1,1));
    startRLSf1 = RLSf1(trialTime < trialStartTime)';
    
    RLSf2 = cell2mat(normalFreqSeries{1,iTrial}(1,2));
    startRLSf2 = RLSf2(trialTime < trialStartTime)';
    
    if ismember(iTrial,onsets(1).(runName))
        F1resolveLow = [F1resolveLow; startRLSf1];
        F2resolveLow = [F2resolveLow; startRLSf2];
    elseif ismember(iTrial,onsets(2).(runName))
        F1resolveHigh = [F1resolveHigh; startRLSf1];
        F2resolveHigh = [F2resolveHigh; startRLSf2];
    end
    
end

loopCount = loopCount+1;
close all
end

%%% Plotting %%%

% close all

%RLS transition data

    
meanF1Lowtrace = mean(F1resolveLow,1);
errorF1Lowtrace = ste(F1resolveLow);

meanF2Lowtrace = mean(F2resolveLow,1);
errorF2Lowtrace = ste(F2resolveLow);

lowStarts(lowStarts > trialStartTime) = []; % remove 'outliers'
meanLowStarts = mean(lowStarts);
errorLowStarts = ste(lowStarts);

meanF1Hightrace = mean(F1resolveHigh,1);
errorF1Hightrace = ste(F1resolveHigh);

meanF2Hightrace = mean(F2resolveHigh,1);
errorF2Hightrace = ste(F2resolveHigh);

meanHighStarts = mean(highStarts);
errorHighStarts = ste(highStarts);

figure
subplot(1,2,1)
title('Resolve Low')
hold on
mseb([(1/512):(1/512):5],[meanF1Lowtrace; meanF2Lowtrace],[errorF1Lowtrace; errorF2Lowtrace],[],1);
vline(meanLowStarts,'k','Start Dominant Percept')
herrorbar(meanLowStarts,0,errorLowStarts,errorLowStarts,'k')
xlabel('Time from trial start (s)');
ylabel('Normalized amplitude');

subplot(1,2,2)
title('Resolve High')
hold on
mseb([(1/512):(1/512):5],[meanF1Hightrace; meanF2Hightrace],[errorF1Hightrace; errorF2Hightrace],[],1);
vline(meanHighStarts,'k','Start Dominant Percept')
herrorbar(meanHighStarts,0,errorHighStarts,errorHighStarts,'k')
xlabel('Time from trial start (s)');
ylabel('Normalized amplitude');

suptitle(parName)
