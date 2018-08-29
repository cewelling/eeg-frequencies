function [cfg_trldef] = ftDefineTrial( EEGfile, numTrials, customDur, event)
% cfg_trldef handles trial definition in the raw EEG data prior to
% pre-processing. Returns trial definition configuration to be used in 
% pre-processing.


% load parameters
analysisParams

if isempty(event)
    event = 'trial';
end

%% Trial Definition
cfg_trldef = [];
cfg_trldef.dataset = EEGfile;
cfg_trldef.trialdef.eventtype = 'STATUS';
cfg_trldef.trialfun = 'ft_trialfun_general';

% Added 1.5 second baseline period after signal for a few specific trials
if strcmp(event, 'trial')
    if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
            || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
            || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
        discard_start = -1.5;
    else
        discard_start = 0;
    end
end

cfg_trldef.trialdef.poststim = trialDur;

if strcmp(event, 'trial')
    cfg_trldef.trialdef.eventvalue = 201:(200+numTrials);
    cfg_trldef.trialdef.poststim = trialDur;
elseif strcmp(event, 'preStim2') || strcmp(event, 'simBase')
    cfg_trldef.trialdef.eventvalue = 201:(200+numTrials);
    cfg_trldef.trialdef.poststim = customDur;
elseif strcmp(event, 'preStim')
    cfg_trldef.trialdef.eventvalue = 45:(44+numTrials);
    cfg_trldef.trialdef.poststim = customDur;
elseif strcmp(event, 'postStim')
    cfg_trldef.trialdef.eventvalue = 99;
    cfg_trldef.trialdef.poststim = customDur;
end

cfg_trldef.trialdef.prestim = discard_start; 

% above: these eventvalues were defined in runScript, which wrote bytes to EEG comp port
% when participant pressed "up" to begin trial

try
    cfg_trldef = ft_definetrial(cfg_trldef);
catch define_trial_error
    cfg_trldef.trialdef.eventvalue
    cfg_trldef.dataset
    rethrow(define_trial_error);
end