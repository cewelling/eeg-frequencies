% Group analysis of binocular rivalry dynamics
% aggregates and analyzes results produced by analysisController.m
% Dependencies: exclude.m

clearvars

%%%%% Switches: Turn on (1) or off (0) %%%%%%%%
groupPlot = 1; % creates group transition plots
fftAnal = 0; % FFT analysis of switch rate
fftCons = 0; % consistency of FFT frequency (trials divided into two groups)
optFFT = 0; % optimize smoothing / filtering for FFT analysis
findCorr = 0; % look for correlations between behavioral and EEG indices
ampClust = 0; % plot group amplitude clusters

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
