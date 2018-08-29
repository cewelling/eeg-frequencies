addpath('./miscFunctions'); %Add misc functions to working directory

%% Clear everything and establish where data is
if exist('E:\Documents\Recorded Data\EEG Feb 2015', 'dir') % location on desktop
    file_directory = 'E:\Documents\Recorded Data\EEG Feb 2015';
elseif exist('D:\Recorded Data', 'dir') % location on laptop
    file_directory = 'D:\Recorded Data';
elseif exist('C:\EEG Data\mit-data', 'dir')
    file_directory = 'C:\EEG Data\mit-data';
elseif exist('/Users/jacksonclee/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results', 'dir')
    file_directory = '/Users/jacksonclee/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results';
elseif exist('/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results', 'dir')
    file_directory = '/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results';
elseif exist('/Users/acspiegel/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results', 'dir')
    file_directory = '/Users/acspiegel/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results';
else
    error('please provide directory where file is stored');
end

plotDir = './plots/';
fftSpectPlotDir = './plots/FFTspect/';
rlsDir = './rlsResults/';
pressDir = '../../runScripts/Results/';

%% Participants
nPars = [1 1]; %# of participants in each group


%% Analysis Variables
discard_start = 0; % how much time should be cut off at beginning
occipitals = [20 26 27 28 29 30 31]; % for 01: [27 28 29 30 31], for 05: [20 26 28 29 30 31] for 08: [27 29 30 31]
electrodes = [1:32]; %[1:32]; %[27, 29, 30,28]; %occipitals

focus_electrode = 29; %27 'LM', 31

noiseThresh = 0.20; %How much signal (prop of max) needs to be in the channel to analyze timepoint
threshInterMod = 0.55; %How much signal (prop of max) needs to be in the channel to analyze timepoint


%% RLS Params

slidingWindowStep = 0.05; %seconds
slidingWindowStep = 0.21; %accomodate ~3 cycles?
slidingWindowWidth = 1; %window specific for each frequency
maxFreq = 30; %hz

%% Diff Score Anal
smoothWindow = 1;


%% Keypress Analysis
%Keypress info
if strcmp(parName, 'stratus01')
    upArrow = 38;
    leftArrow = 37;
    rightArrow = 39;
else
    upArrow = 104;
    leftArrow = 100;
    rightArrow = 102;
end

pressLatency = 0.4; %0.2

minEpoch = 0.5; %minimum cutoff for the duration of perceptual epochs

