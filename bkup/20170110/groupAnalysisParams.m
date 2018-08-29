% Parameters for group analysis 

%% Directories-------------------------------------------------------------

% location of behavioral results
bResultsDir = '/Users/acspiegel/Dropbox (MIT)/Projects/EEG_Frequencies/analScripts/Behavior/analResults/';

% Location of Analysis Results
indicesDir = ['./indices/'];
corrPlotDir = './plots/correlations/';
groupPlotDir = './plots/Group/';

%% Group participants------------------------------------------------------

groupParNums = {101:138, 01:20};
groupCodes = {'stratus', 'cumulus'};
analGroupIDs = {'stratus', 'cumulus'};

% Create label indicating participants analyzed
label = [];
for iGroup = 1:length(groupParNums)
    groupLabel = [groupCodes{iGroup} num2str(groupParNums{iGroup}(1)) '-' num2str(groupParNums{iGroup}(end))];
    if iGroup == 1
        label = [label groupLabel];
    else
        label = [label ', ' groupLabel];
    end
end

%% Types of runs-----------------------------------------------------------

runTypes = {'dartRival','dartSim','marzRival'}; %'dartSim', 'marzRival'
runIndices = [1:3]; % generally 1:3

%% Gating parameters ------------------------------------------------------
% which participants get included in the analysis?

minSNR = 1;
%minRedEpochs = ;
%minGreenEpochs = ;
%number of runs per participant?
% Produce lists of included and excluded participants to be
% included in the folder for this analysis

%% SNR parameters ---------------------------------------------------------

occipitals = struct('nums', [20 26 27 28 29 30 31], 'name', 'occipitals');
focusElectrode = struct('nums', 29, 'name', 'Oz');
allElectrodes = struct('nums', 1:32, 'name', 'allElectrodes');

% Choose set of electrodes to analyze:
electrodeSet = occipitals; 

noiseWindowHalf = 3; % in freq bins for SNR calculation--1 Hz = 3 bins
FOI = 8.5; % hz 5.666 or 8.5

%% Statistics parameters---------------------------------------------------

% Which types of rivalry runs to analyze
dartboards = 1:3; % numbers associated with darboard runs in legends
marz = 4:6;
stimType = marz;

%% Plotting parameters-----------------------------------------------------
transHalf = 5; %+/- from the time of the button press %1.5
peakHalf = 5; % in seconds, +/- around peak
sampRate = 512;
minInstances = 2;

%% Transition parameters---------------------------------------------------
allTrans = 0; % include transitions that go through mixed states in low to high and high to low?
buttonPress = transHalf* sampRate + 1; % index of button press in transitions
