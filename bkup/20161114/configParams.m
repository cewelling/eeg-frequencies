function [cfg_trldef, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,electrodes)

% Trial Definition
cfg_trldef.dataset = fileID;
cfg_trldef.trialdef.eventtype = 'STATUS';
cfg_trldef.trialfun = 'ft_trialfun_general';
cfg_trldef.trialdef.prestim = discard_start;
cfg_trldef.trialdef.poststim = trial_dur;
cfg_trldef.trialdef.eventvalue = 201:(200+trialNum); %216 JF 
cfg_trldef = ft_definetrial(cfg_trldef);

% Pre-processing
cfg_preproc = cfg_trldef;
cfg_preproc.channel = 1:32; %1:32

cfg_preproc.continuous = 'yes';
cfg_preproc.demean = 'yes';
cfg_preproc.detrend = 'yes';
cfg_preproc.reref = 'yes';
% cfg_preproc.implicitref = 'RM';
cfg_preproc.refchannel = 'all'; %{'01', 'Oz', 'Iz', 'POz'}; %%electrodes

cfg_preproc.bpfilter = 'yes'; % bandpass filter?
cfg_preproc.bpfreq = [2 70]; % bandpass frequencies 
%cfg_preproc.bpfreq = [1 30]; % bandpass frequencies 

cfg_preproc.bsfilter = 'yes'; % bandstop filter?
cfg_preproc.bsfreq = [59 61]; % bandstop frequencies (US = 60)

% cfg_preproc.hpfilter = 'yes';
% cfg_preproc.hpfreq = 2;