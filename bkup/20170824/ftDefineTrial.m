function [cfg_trldef] = ftDefineTrial( EEGfile, numTrials, trialDur)
% cfg_trldef defines segments (trials) in the raw EEG data to be analyzed.
% Returns trial definition configuration to be used in pre-processing.

% General trial definition parameters
cfg_trldef = [];
cfg_trldef.dataset = EEGfile;
cfg_trldef.trialdef.eventtype = 'STATUS';
cfg_trldef.trialfun = 'ft_trialfun_general';

% Added 1.5 second baseline period after signal for a few specific trials
if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
    discard_start = -1.5;
else
    discard_start = 0;
end

% eventvalues were defined in runScript, which wrote bytes to EEG comp port
% when participant pressed "up" to begin trial (the trigger)
cfg_trldef.trialdef.eventvalue = 201:(200+numTrials);

% beginning of trial (relative to trigger)
cfg_trldef.trialdef.prestim = discard_start;
% end of trial (relative to start of trial)
cfg_trldef.trialdef.poststim = trialDur;

% Run trial definition
try
    cfg_trldef = ft_definetrial(cfg_trldef);
catch define_trial_error
    cfg_trldef.trialdef.eventvalue
    cfg_trldef.dataset
    rethrow(define_trial_error);
end

end