function [ all_data ] = ftConfigParams( parName, runName, EEGfile, numTrials, trialDur, date )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
analysisParams

%% Trial Definition
cfg_trldef = [];
cfg_trldef.dataset = EEGfile;
cfg_trldef.trialdef.eventtype = 'STATUS';
cfg_trldef.trialfun = 'ft_trialfun_general';

% Added 1.5 second baseline period after signal for a few specific trials
if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
    discard_start = 1.5;
end

cfg_trldef.trialdef.prestim = discard_start + 1; %discard_start;
cfg_trldef.trialdef.poststim = trialDur;
cfg_trldef.trialdef.eventvalue = [201:(200+numTrials)];
% these eventvalues were defined in runScript, which wrote bytes to EEG comp port
% when participant pressed "up" to begin trial

% Mistakes-----------------------------------------------------------------

% Did not save first two trials
if ~isempty(strfind(EEGfile, 'cumulus05_marzRival1'))
    cfg_trldef.trialdef.eventvalue = 203:206;
end

%--------------------------------------------------------------------------

try
    cfg_trldef = ft_definetrial(cfg_trldef);
catch define_trial_error
    cfg_trldef.trialdef.eventvalue
    cfg_trldef.dataset
    rethrow(define_trial_error);
end

%% Pre-pre-processing (Identify noisy electrodes)
if ~exist(['badElecs/' parName '_' runName '.mat'], 'file')
    cfg_preproc = cfg_trldef;
    cfg_preproc.channel = 'all';
    cfg_preproc.continuous = 'yes';
    cfg_preproc.demean    = 'yes'; % perform baseline correction
    
    % for a few trials, baseline window comes after signal
    if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
            || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
            || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
        cfg_preproc.baselinewindow = [0 1.5];
    end
    
    % Started adding a baseline window before signal
    date = datetime(date, 'Format','yyyy-MM-dd');
    if date > datetime('2016-12-02','Format','yyyy-MM-dd');
        %cfg_preproc.baselinewindow = [-1.5 0];
    end
    
    cfg_preproc.detrend = 'yes'; % remove linear trend from the data
    
    cfg_preproc.reref = 'no';
    
    cfg_preproc.hpfilter = 'yes'; % highpass filter
    cfg_preproc.hpfreq = 2; % highpass frequency in hertz
    
    checkNoise_data = ft_preprocessing(cfg_preproc);
    
    %% Find noisy electrodes
    
    cfg_rejectvis = [];
    cfg_rejectvis.method = 'summary';
    cfg_rejectvis.keepchannel = 'nan';
    pruned_data = ft_rejectvisual(cfg_rejectvis,checkNoise_data);
    
    % store noisy electrodes
    badElecs = [];
    trialData = pruned_data.trial{1}; % all trials have the same electrodes eliminated, only need to check one
    for iElec = 1:size(trialData,1);
        if isnan(trialData(iElec, 1))
            badElecs = [badElecs iElec];
        end
    end
    save(['badElecs/' parName '_' runName], 'badElecs');
    
else
    load(['badElecs/' parName '_' runName '.mat']);
end

%% Preprocessing

cfg_preproc = cfg_trldef;
cfg_preproc.channel = 'all';
cfg_preproc.continuous = 'yes';
cfg_preproc.demean    = 'yes'; % perform baseline correction

% for a few trials, baseline window comes after signal
if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
    cfg_preproc.baselinewindow = [0 1.5];
end

% Started adding a baseline window before signal
date = datetime(date, 'Format','yyyy-mm-dd');
if date > datetime('2016-12-02','Format','yyyy-mm-dd');
    %cfg_preproc.baselinewindow = [-1.5 0];
end

cfg_preproc.detrend = 'yes'; % remove linear trend from the data

% re-reference
cfg_preproc.reref = 'yes';
refElecs = num2cell([1:32]'); %1:32
for i = 1:length(refElecs)
    refElecs{i} = num2str(refElecs{i});
end
refElecs = [refElecs; 'veog'; 'heog'; 'lm'; 'rm']; %, 'veog'; 'heog'; 'lm'; 'rm'}; % veog and heog are mastoid channels, lm and rm are eyes

% preclude noisy electrodes from being included as a reference
goodIndices = 1:length(refElecs);
for i = 1:length(badElecs)
    if badElecs(i) < 33
        badElec = num2str(badElecs(i));
    elseif badElecs(i) == 65
        badElec = 'veog';
    elseif badElecs(i) == 66
        badElec = 'heog';
    elseif badElecs(i) == 67
        badElec = 'lm';
    elseif badElecs(i) == 68
        badElec = 'rm';
    else
        badElec = [];
    end
    badIndex = find(strcmp(refElecs, badElec));
    if ~isempty(badIndex)
        goodIndices = goodIndices(find(goodIndices ~= badIndex));
    end
end
refElecs = refElecs(goodIndices); 

 cfg_preproc.refchannel = refElecs; % common average or median reference
% cfg_preproc.refchannel = 'all'; % common average reference
 cfg_preproc.refmethod = 'median'; % median %cer 11/16
% cfg_preproc.implicitref = 'M1';            % the implicit (non-recorded) reference channel is added to the data representation
% cfg_preproc.refchannel     = {'M1', 'M2'}; % the average of these channels is used as the new reference, note that channel '53' corresponds to the right mastoid (M2)

cfg_preproc.hpfilter = 'yes'; % highpass filter
cfg_preproc.hpfreq = 2; % highpass frequency in hertz

all_data = ft_preprocessing(cfg_preproc);

%% Remove trials that aren't real 
% A mistake in our runScripts was sending some misplaced trial-marking signals to
% the EEG computer 

% stratus138
if ~isempty(strfind(EEGfile, 'stratus138_dartRival1')) || ~isempty(strfind(EEGfile, 'stratus138_dartSim3')) ...
        || ~isempty(strfind(EEGfile, 'stratus138_marzRival1'))
    all_data = removeTrials(all_data, numTrials, 3);
elseif ~isempty(strfind(EEGfile, 'stratus138_dart')) || ~isempty(strfind(EEGfile, 'stratus138_marzRival2'))
    all_data = removeTrials(all_data, numTrials, 2, 4);
elseif ~isempty(strfind(EEGfile, 'stratus138'))
    all_data = removeTrials(all_data, numTrials, 2);
end

% cumulus15
if ~isempty(strfind(EEGfile, 'cumulus15_dartSim1')) || ~isempty(strfind(EEGfile, 'cumulus15_dartSim2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus15_marzRival3'))
    all_data = removeTrials(all_data, numTrials, 3);
elseif ~isempty(strfind(EEGfile, 'cumulus15_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus15_marzRival'))
    all_data = removeTrials(all_data, numTrials, 2, 4);
elseif ~isempty(strfind(EEGfile, 'cumulus15'))
    all_data = removeTrials(all_data, numTrials, 2);
end

% cumulus05
if ~isempty(strfind(EEGfile, 'cumulus05_dartRival2')) || ~isempty(strfind(EEGfile, 'cumulus05_dartSim2'))
    all_data = removeTrials(all_data, numTrials, 2, 4);
elseif ~isempty(strfind(EEGfile, 'cumulus05_dartSim3')) || ~isempty(strfind(EEGfile, 'cumulus05_marzRival2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus05_marzRival3'))
    all_data = removeTrials(all_data, numTrials, 2);
end

% cumulus06
if ~isempty(strfind(EEGfile, 'cumulus06_marzRival1')) 
    all_data = removeTrials(all_data, numTrials, 3);
elseif ~isempty(strfind(EEGfile, 'cumulus06_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus06_dartSim1')) ...
        || ~isempty(strfind(EEGfile, 'cumulus06_dartSim3')) || ~isempty(strfind(EEGfile, 'cumulus06_marzRival'))
    all_data = removeTrials(all_data, numTrials, 2, 4);
elseif ~isempty(strfind(EEGfile, 'cumulus06'))
    all_data = removeTrials(all_data, numTrials, 2);
end


end

