% Group analysis of binocular rivalry dynamics
% aggregates and analyzes results produced by analysisController.m
% Dependencies: exclude.m

clearvars

%%%%% Switches: Turn on (1) or off (0) %%%%%%%%
groupPlot = 0; % creates group transition plots
fftAnal = 0; % FFT analysis of switch rate
fftCons = 0; % consistency of FFT frequency (trials divided into two groups)
optFFT = 0; % optimize smoothing / filtering for FFT analysis
findCorr = 1; % look for correlations between behavioral and EEG indices
ampClust = 0; % plot group amplitude clusters

%%%%% Machine learning switches: Turn on (1) or off (0) %%%%%%%%%%%%%%%%%%%

% edit parameters in svmParams.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classifyPercept = 0; % general classification script (percepts, not diag)
% runSVMClassifier.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% edit parameters in individual files:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crossVal = 0; % Cross-validate by training on 3 excluded participants, testing on 1
% SVMcrossVal.m
airplane = 0; % Optimize testing segment and latency
% runSVMAirplaneMatrices.m
diagnose = 0; % Classify autists vs. controls
% SVMclassAvsC.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prep
%groupNum = 1
% Exclude participants (see exclude.m for list of participants excluded)
exclude

% load parameters
groupAnalysisParams % must be loaded AFTER running exclude.m

%% Group plotting

if groupPlot == 1
    makeGroupPlot;
end

%% Group FFT oscillation analysis

% Compare across groups
if fftAnal == 1
    avgFFT;
end

% Test-retest reliability
if fftCons == 1
    FFTconsistency;
end

% Optimize smoothing
if optFFT == 1
    optSmoothing;
end

%% Statistics

if findCorr == 1
    indexAnal2;
end

%% Group amplitude clustering

if ampClust == 1
    groupAmpClusters;
end

%% Machine learning

if classifyPercept 
    runSVMclassifyPercept;
end

if crossVal 
    runSVMcrossVal;
end

if airplane 
    runSVMAirplaneMatrices;
end

if diagnose 
    classAvsC;
end
