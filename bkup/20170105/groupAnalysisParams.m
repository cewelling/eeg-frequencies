% Parameters for group analysis 

%% Directories-------------------------------------------------------------

% location of behavioral results
bResultsDir = '/Users/acspiegel/Dropbox (MIT)/Projects/EEG_Frequencies/analScripts/Behavior/analResults/';

% Location of Analysis Results
indicesDir = ['./indices/'];

% Legend location
legendDir = ['./pilotlegends/'];

%% Group participants------------------------------------------------------

groupParNums = {101:137, 03:07};
groupCodes = {'stratus', 'cumulus'};
analGroupIDs = {'stratus', 'cumulus'};

%% Types of runs-----------------------------------------------------------

runIndices = [7 8 9]; % Corresponds to legend; typical cases below
% [1 2 3] for dartboard rivalry
% [4 5 6] for dartboard sim
% [7 8 9] for marz rivalry

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
sampRate = 512;
minInstances = 20;

%% Transition parameters---------------------------------------------------
buttonPress = transHalf* sampRate + 1; % index of button press in transitions
