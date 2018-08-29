% -------------------------------------------------------------------------
% analysisParams.m
% description: parameters for SSVEP binocular rivalry anlysis
% Use with: analysisController.m, runFFT.m
% -------------------------------------------------------------------------

%% Plot settings
FFTplotOrNot = 'no'; 
RLSplotOrNot = 'no'; 
transPlotOrNot = 'yes'; % transitions for each run
PkPlotOrNot = 'no'; % peaks overlayed on RLS spectra
avgPksPlotOrNot = 'no'; % peaks for each run
parPlotOrNot = 'no';
clustPlotOrNot = 'no';

%% Directories-------------------------------------------------------------

% Location of Raw Data
if exist('/Users/acspiegel/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/rawData', 'dir')
    eegDir = '/Users/acspiegel/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/rawData/EEG/'; 
    keyPressDir = '/Users/acspiegel/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/rawData/Behavior/'; 
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
blineDir = './baselineTest/';

% Date format in file names
dateFormat = 'yyyy-mm-dd';

%% Experiment parameters---------------------------------------------------
group = 0;
if group == 0
    parNums = [105:140]; %[105:140]; % ex/ Stratus100's parNum is 100
    group = 'stratus'; % 'stratus' or 'cumulus'
else
    parNums = [1:29]; %[1:29]; % ex/ Stratus100's parNum is 100
    group = 'cumulus'; % 'stratus' or 'cumulus'
end
runTypes = {'dartRival'}; %'dartSim', 'marzRival'
runIndices = [1:3]; % generally 1:3

stimFreqs = [5.666 8.5]; % hz
harFreqs = [2*stimFreqs(1) 2*stimFreqs(2)]; % harmonics
imFreqs = [stimFreqs(2) - stimFreqs(1) stimFreqs(1) + stimFreqs(2)]; % intermodulation frequencies

analFreqs = stimFreqs; % stimFreqs;

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

% Choose sets of electrodes to analyze:
electrodeSets = [focusElectrode]; %focusElectrode; %[occipitals]; 
numElecs =  1; %1; %7;

discard_start = 0; %0; % time to cut off of beginning of EEG data
noiseWindowHalves = [6]; % in freq bins for SNR calculation--1 Hz = 3 bins 

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
electrodeSet = allElectrodes; % Choose set of electrodes to analyze:
smoothingPts = 150; % __ pt moving average to smooth RLS amplitudes

% gaussian kernel
smoothingWin = 1000; %150 
smoothingStd = 500; %20

% haoran's smoothing scripts (modified to accomodate just a vector of data,
% rather than any 2 dimensional matrix)
smooth_ksize = 1225; % must be an odd integer 
smooth_sd = 1; %std is expressed as a fraction of (smooth_ksize + 1) /2 

% Transition Analysis------------------------------------------------------

% Time requirements for consideration in transition analysis
dom1Min = 0.5; %1.25;
dom1Max = inf;
dom2Min = 0.5;
dom2Max = inf;
mixedMin = .75; %.2 %1;  
mixedMax = 1.3; %.75
spaceMin = 0;
spaceMax = inf;
noMixed = 0;
gapMax = 0.25; % max time allowed between button-reported states

transHalf = 5; %+/- from the time of the button press %1.5
buttonPress = transHalf*512 + 1; % index of button press in transitions

normType = 'mean'; %'noStimBase'; %'freqWin'; % or 'z' or 'none' or 'norm' or 'mean'
blineLoc = 'trial'; % or 'preStim' or 'postStim' or 'trial' (other freqs during the trial)

% for modulation index calculation
minSNR = 3;
minInstMod = 5; % min number of transitions for each frequency, run type, and transition type

% Peak Analysis------------------------------------------------------------

% pick peaks without relying on perceptual report?
blindPeak = 0;
maxPwidth = 5;
minPwidth = 1;

% alternate (simplistic) peak analysis?
altPeak = 0;

peakHalf = 5; % in seconds, +/- around peak
minPkInstances = 5;

% Participant Plotting-----------------------------------------------------

minInstances = 5;
commonElecs = 1;
