function [ preproc_data ] = ftPreProc( cfg_trldef, parName, EEGfile, date, reref, lowBound, highBound,cfg_preproc)
% ftPreProc returns pre-processed data with or without referencing, as 
% indicated by the argument reref
% lowBound and highBound are the bounds of the band pass filter

cfg_preproc = cfg_trldef;


%% Baseline Window: Not currently using this, but good to know which participants have it
% for a few trials, baseline window comes after signal
% if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
%         || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
%         || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
%     %cfg_preproc.baselinewindow = [0 1.5];
% end

% Started adding a baseline window before signal
% date = datetime(date, 'Format','yyyy-MM-dd');
% if date > datetime('2016-12-02','Format','yyyy-MM-dd');
%     %cfg_preproc.baselinewindow = [-1.5 0];
% end


%% Re-referencing
if reref
    cfg_preproc.reref = 'yes';
    refElecs = 1:32;
    cfg_preproc.refchannel = refElecs; % common average or median reference
    cfg_preproc.refmethod = 'avg';
else
    cfg_preproc.reref = 'no';
end

%% Filtering
cfg_preproc.bpfreq = [lowBound highBound]; %[4 10]; %[2 70]; % bandpass frequencies

%% Run pre-processing
preproc_data = ft_preprocessing(cfg_preproc);


end

