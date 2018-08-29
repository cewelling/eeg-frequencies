% -------------------------------------------------------------------------
% analysisParams.m
% description: parameters for SSVEP binocular rivalry anlysis
% Use with: analysisController.m, runFFT.m
% -------------------------------------------------------------------------

%% Plot settings
FFTplotOrNot = 'no'; 
RLSplotOrNot = 'no'; 
transPlotOrNot = 'no'; % transitions for each run
PkPlotOrNot = 'no'; % peaks overlayed on RLS spectra
avgPksPlotOrNot = 'no'; % peaks for each run
parPlotOrNot = 'no';
clustPlotOrNot = 'no'; % amplitude clustering segment plots
segPlotOrNot = 'no'; % plots of machine learning trials

% plotting properties
setFigProps

color3 = [0.15 0.15 0.15];
colorA = [110 22 48]/255;
colorC = [10 77 97]/255;

%% Directories-------------------------------------------------------------

% Location of Raw Data
if exist('../../runScripts/rawData', 'dir')
    eegDir = '../../runScripts/rawData/EEG/'; 
    keyPressDir = '../../runScripts/rawData/Behavior/'; 
elseif exist('/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/rawData', 'dir')
    eegDir = '/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/rawData/EEG/'; 
    keyPressDir = '/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/rawData/EEG/';   
else
    error('please provide directory where raw data is stored');
end

% Location of Analysis Results
indicesDir = ['./indices/'];
plotDir = './plots/';
fftSpectPlotDir = './plots/FFTspect/';
transPlotDir = './plots/transitions/';
blineDir = '.pre-processing/baselineTest/';

% Date format in file names
dateFormat = 'yyyy-mm-dd';

%% Experiment parameters---------------------------------------------------
group = 0; %0 controls, 1 for asc
if group == 0
    parNums = [106 112 131 132]; %[105:121 123:128 130:132 134 136:140]; %[105:140]; % ex/ Stratus100's parNum is 100
    group = 'stratus'; % 'stratus' or 'cumulus'
else
    parNums = [1:22 24:28]; %[1:28]; % ex/ Stratus100's parNum is 100
    group = 'cumulus'; % 'stratus' or 'cumulus'
end
runTypes = {'dartRival' 'dartSim'}; %'dartSim', 'marzRival'
runIndices = [1:3]; % generally 1:3

stimFreqs = [5.666 8.5]; % hz
harFreqs = [2*stimFreqs(1) 2*stimFreqs(2)]; % harmonics
imFreqs = [2*stimFreqs(2) - stimFreqs(1) stimFreqs(1) + stimFreqs(2)]; % intermodulation frequencies

analFreqs = imFreqs; % stimFreqs; % harFreqs; % imFreqs

if sum(analFreqs == stimFreqs) == 2
    analFreqLabel = 'principles';
elseif sum(analFreqs == harFreqs) == 2
    analFreqLabel = 'harmonics';
elseif sum(analFreqs == imFreqs) == 2
    analFreqLabel = 'im'; 
end

numTrials = 6;
trialDur = 30; % in seconds

sampRate = 512; % EEG sample rate (per second)

%% Analysis variables------------------------------------------------------

% FFT----------------------------------------------------------------------
% Note that there are a few settings that must be changed within runFFT
% (ex/ Hanning taper vs. multitaper, amount of padding, etc.)

occipitals = struct('nums', [20 26 27 28 29 30 31], 'name', 'occipitals');
focusElectrode = struct('nums', 29, 'name', 'Oz');
focusElectrodes = struct('nums', [20 29], 'name', 'Oz and POz');
allElectrodes = struct('nums', 1:32, 'name', 'allElectrodes');
mastoids = struct('nums', 65, 'name', 'mastoids');
custom = struct('nums', [1:10], 'name', '1-10');

% Choose electrodes to analyze:
electrodeSet = occipitals; %focusElectrode; %occipitals;
snrElecs = focusElectrode.nums; % electrodes used for SNR exclusion purposes
%avgOrInd = 'ind'; % average electrodes? or analyze each one individually? 

discard_start = 0; % time to cut off of beginning of EEG data
noiseWindowHalves = [15]; % in freq bins for SNR calculation--1 Hz = 30 bins 

% Keypress Analysis--------------------------------------------------------

%Keypress info
upArrow = 104;
leftArrow = 100;
rightArrow = 102;

%For simulation schedule
leftCue = 1;
rightCue = 2;

minEpochDur = 0.5; %0.5; % in seconds

% RLS Analysis-------------------------------------------------------------

cutoffWindow = 1000; % in ms
smoothingPts = 150; % __ pt moving average to smooth RLS amplitudes

% haoran's smoothing scripts (modified to accomodate just a vector of data,
% rather than any 2 dimensional matrix)
smooth_ksize = 1225; % must be an odd integer 
smooth_sd = 1; %std is expressed as a fraction of (smooth_ksize + 1) /2 

% Transition Analysis------------------------------------------------------

% use sim schedule for simulation trials?
useSimSchedule = 1;

% Time requirements for consideration in transition analysis
dom1Min = 0.5; %0.5 %1.25; %cutoff for shortest dominant state to be included
dom1Max = inf;
dom2Min = 0.5; %0.5 %cutoff for shortest dominant state to be included
dom2Max = inf;
mixedMin = 0; %.2 %1;  
mixedMax = inf; %.75
noMixed = 0; % if this is 1, do not include transitions with a mixed component
gapMax = 0.25; % max time allowed between button-reported states

transHalf = 5; %+/- from the time of the button press %1.5
buttonPress = transHalf*512 + 1; % index of button press in transitions

normType = 'mean'; %'noStimBase'; %'freqWin'; % or 'z' or 'none' or 'norm' or 'mean'
blineLoc = 'trial'; % or 'preStim' or 'postStim' or 'trial' (other freqs during the trial)

% for modulation index calculation
minSNR = 2;
minInstMod = 5; % min number of transitions for each frequency, run type, and transition type

% Peak Analysis------------------------------------------------------------

% pick peaks without relying on perceptual report?
blindPeak = 0;
maxPwidth = 5;
minPwidth = 1;

% alternate (simplistic) peak analysis? Just choose the midpoint between
% start and end of dominant state (according to perceptual report)
altPeak = 0;

peakHalf = 5; % in seconds, +/- around peak
minPkInstances = 5;

% Participant Plotting-----------------------------------------------------

minInstances = 5;

% AmpClusters Analysis-----------------------------------------------------

latency = round(.8*sampRate); % predicted latency in data points
segFrac = 1; % proportion of epoch used for amplitude clustering

%% Machine Learning Features

amp = 1; % amplitude
d1dt = 1; % 1st derivative
d2dt = 1; % 2nd derivative

% frequencies
freqTypes = {'principles' 'harmonics' 'im'}; % 'harmonics' 'im'}; %{'principles' 'harmonics' 'im'};

% channels
channels = occipitals.nums;

% transition parameters
svmTransHalf = 3;
