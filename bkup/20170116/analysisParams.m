% -------------------------------------------------------------------------
% analysisParams.m
% description: parameters for SSVEP binocular rivalry anlysis
% Use with: analysisController.m, runFFT.m
% -------------------------------------------------------------------------

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

% Date format in file names
dateFormat = 'yyyy-mm-dd';

%% Experiment parameters---------------------------------------------------

parNums = [125]; % ex/ Stratus100's parNum is 100
group = 'stratus'; % 'stratus' or 'cumulus'
runTypes = {'dartSim'}; %'dartSim', 'marzRival'
runIndices = [2:3]; % generally 1:3

stimFreqs = [5.666 8.5];
%imFreqs = linear combos of stimFreqs1 & 2; %intermodulation frequencies
numTrials = 6;
trialDur = 30; % in seconds

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
electrodeSets = [allElectrodes]; 
numElecs = 20;

discard_start = 0; %0; % time to cut off of beginning of EEG data
noiseWindowHalves = [6]; % in freq bins for SNR calculation--1 Hz = 3 bins 

% Create plots?
FFTplotOrNot = 'no'; % 'yes' or 'no'

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
smooth_sd = 1; %std is expressed as a fraction of smooth_ksize + 1 /2 

% Create plots?
RLSplotOrNot = 'no'; % 'yes' or 'no'

% Transition Analysis------------------------------------------------------

% Include transitions that go through mixed states?
allTrans = 1;

% Time requirements for consideration in transition analysis
domMin = 1.25;
mixedMin = 1;
gapMax = 0.25; % max time allowed between button-reported states

transHalf = 5; %+/- from the time of the button press %1.5
buttonPress = transHalf*512 + 1; % index of button press in transitions
normType = 'mean'; % or 'z' or 'none' or 'norm'
% Create plots for individual runs?
transPlotOrNot = 'no'; % 'yes' or 'no'

% for modulation index calculation
minSNR = 3;
minInstMod = 5; % min number of transitions for each frequency, run type, and transition type

% Peak Analysis------------------------------------------------------------

% alternate (simplistic) peak analysis?
altPeak = 0;

PkPlotOrNot = 'yes';
peakHalf = 5; % in seconds, +/- around peak
avgPksPlotOrNot = 'yes';
minPkInstances = 5;

% Participant Plotting-----------------------------------------------------

minInstances = 5;
commonElecs = 2;
