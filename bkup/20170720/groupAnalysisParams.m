% Parameters for group analysis 

%% Directories-------------------------------------------------------------

% location of raw and behavioral data
if exist('../../runScripts/rawData', 'dir')
    eegDir = '../../runScripts/rawData/EEG/'; %eeg data
    keyPressDir = '../../runScripts/rawData/Behavior/'; %raw behavioral data
    bResultsDir = '../../analScripts/Behavior/analResults/'; %behavioral data
else
    error('please provide directory where raw data is stored');
end

% Date format in file names
dateFormat = 'yyyy-mm-dd';

% Location of Analysis Results
indicesDir = './indices/';
corrPlotDir = './plots/correlations/';
barPlotDir = './plots/barGraphs/';
groupTransPlotDir = './plots/transitions/Group/';

%% Group participants------------------------------------------------------

groupParNums = {[105:138], [1:28]}; %{[105:138], [1:28]}; %{[105:138], [1:28]}; %{[106 112 131 132], []}; %{[105:124 126:138], [1:29]};%{[105:138], [01:13 15:20]}; %, 01:20};
groupCodes = {'stratus','cumulus'}; %{'stratus','cumulus'};
analGroupIDs = {'stratus', 'cumulus'}; %{'stratus', 'cumulus'};

% exclude participants
load('pre-processing/stratExcluded.mat')
load('pre-processing/cumExcluded.mat')
groupParNums{1} = setdiff(groupParNums{1},stratExcluded);
groupParNums{2} = setdiff(groupParNums{2},cumExcluded);

% Create label indicating participants analyzed
label = {};
% for iGroup = 1:length(groupParNums)
%     groupLabel = [groupCodes{iGroup} num2str(groupParNums{iGroup}(1)) '-' num2str(groupParNums{iGroup}(end))];
%     label{iGroup} = groupLabel;
% end

%% Types of runs-----------------------------------------------------------

runTypes = {'dartRival','dartSim'}; %'dartSim', 'marzRival'
runIndices = [1:3]; % generally 1:3

%% Gating parameters ------------------------------------------------------

% which participants get included in the analysis?
%minSNR = 3;

%% Electrodes -------------------------------------------------------------

% Define electrode sets
occipitals = struct('nums', [20 26 27 28 29 30 31], 'name', 'occipitals');
focusElectrode = struct('nums', 29, 'name', 'Oz');
focusElectrodes = struct('nums', [20 29], 'name', 'Oz and POz');
allElectrodes = struct('nums', 1:32, 'name', 'allElectrodes');
mastoids = struct('nums', 65, 'name', 'mastoids');
custom = struct('nums', [1:10], 'name', '1-10');

%% Frequencies ------------------------------------------------------------

% Define frequencies
stimFreqs = [5.666 8.5]; % hz
harFreqs = [2*stimFreqs(1) 2*stimFreqs(2)]; % harmonics
imFreqs = [2*stimFreqs(2) - stimFreqs(1) stimFreqs(1) + stimFreqs(2)]; % intermodulation frequencies

% Define frequency labels
stimFreqLabel = 'principles';
harFreqLabel = 'harmonics';
analFreqLabel = 'im';

%% SNR parameters----------------------------------------------------------

%noiseWindowHalf = 3; % in freq bins for SNR calculation--1 Hz = 3 bins
%FOI = 8.5; % hz 5.666 or 8.5

%% Group plot parameters---------------------------------------------------
sampRate = 512;
minInstances = 4;
minPars = 5;
normType = 'mean'; %'none'; % 'freqWin'

% Choose electrodes to analyze:
electrodeSet = focusElectrode; %focusElectrode; %occipitals;

% Choose frequencies to analyze:
analFreqs = stimFreqs; % stimFreqs; % harFreqs; % imFreqs
analFreqLabel = stimFreqLabel;

% Transition parameters---------------------------------------------------
transHalf = 5; %+/- from the time of the button press %1.5
allTrans = 0; % include transitions that go through mixed states?
buttonPress = transHalf* sampRate + 1; % index of button press in transitions

% use sim schedule
useSimSchedule = 0;

% Peak parameters---------------------------------------------------------
peakHalf = 5; % in seconds, +/- around peak
altPeak = 0; % if 1, simplistic peak analysis (line up middle of epoch)
peakI = peakHalf* sampRate + 1;

% Alignment of normalized group plot
alignToCross = 0; % otherwise, align to button press

%% FFT analysis parameters-------------------------------------------------

%% Correlation parameters--------------------------------------------------

%% Plotting settings-------------------------------------------------------
setFigProps
color3 = [0.15 0.15 0.15];
colorA = [110 22 48]/255;
colorC = [10 77 97]/255;
