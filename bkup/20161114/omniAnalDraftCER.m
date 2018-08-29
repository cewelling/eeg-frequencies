%%Meta analysis script for EEG binocular rivalry data
%%Calls a sequence of analysis steps, which can be turned on and off under "switches"
%Depends on Parameters stored in defaultAnalParams.m
%Also depends on RLS parameters stored in configParams.m
%March 2016

clear all; close all; clc;
% setFigProps

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Set up: Things to change %%%%%%%%%%%%%%%
analRunNum = 1;
parName = 'stratus105';
% date = '4_15';
% fileNames = dir(['../../runScripts/Results/' parName '*' date '*']);
fileNames = dir(['../../runScripts/Results/' parName '*']);
legend = importdata(['pilotlegends/' parName '.txt']); %pilot-' date '-caroline.txt
stimulation_frequencies = [14.16 17]; %[17 21.25] Note: this is redefined later....
intermodulation_frequencies = 11.3;

%%%%% Switches: Turn on (1) or off (2) %%%%%%%%
fftAnal = 1; %plain fft
rlsAnal = 0; %rls analysis
diffAnal = 0; %difference scores between two EEG freqs
keyAnal = 0; %parses keypress data
alignAnal = 0; %aligns keypress with EEG

%%%%% Set up: Always done %%%%%%%%%%%%%%%%%%%%%
fileNames = {fileNames.name}; filenames{1}(1).name = char(fileNames(analRunNum));
stimulation_frequencies = legend.data(analRunNum,1:2);
numTrials = legend.data(analRunNum,3);
trialDur = legend.data(analRunNum,4); 
defaultAnalParams %Load default anal params

% stimulation_frequencies = [7.0833 8.5000];

stimulation_frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Analysis             %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run FFT
if fftAnal == 1
    [maxResponseElecs, maxIndices] = runFFT(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend)
end

%% Run RLS
if rlsAnal == 1
    runRLS(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend, maxIndices)
end

%% Compute Frequency Difference Scores
if diffAnal == 1
    runDiffAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName)
end

%% Parse Keypress Data
if keyAnal == 1
    omniData = runKeyAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,minEpoch)
    
    pressDurations(omniData,parName,analRunNum)
    
    plotRLSspectra_meanAcrossTrials(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend,omniData)

end

%% Compute Keypress Alignment with EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alignAnal == 1
    omniData = runKeyAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,minEpoch)
    
    crossCorDiffPress(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,omniData)
end

% [simLatency, latencyError] = determineLatency(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName)

%% Working

omniData = runKeyAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,minEpoch);

% [slowDur, fastDur, mixedDur] = pressDurations(omniData,parName,analRunNum);
[slowDur, fastDur] = pressDurations(omniData,parName,analRunNum);

compareRLS(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend, maxIndices, omniData)
% [multiFreqSeries, trialTime] = compareRLSfast(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend, omniData);

[transitions] = compareEpochs(omniData, multiFreqSeries, trialTime);

% [maxResponseElecs, maxIndices] = runFFT(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend)
% 
% runRLS(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend, maxIndices)

