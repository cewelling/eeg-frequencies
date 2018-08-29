function clean_data = removeArtifacts(data,cfg_trldef,interactive)

cfg = [];

cfg.trl = cfg_trldef.trl; %cer .trl

cfg.artfctdef.reject = 'nan'; % 'partial' to remove sections of trial
cfg.artfctdef.minaccepttim = 5; %minimum lenth of data to retain (seconds)

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = 'veog'; %use the external eye electrode
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

% feedback
cfg.artfctdef.zvalue.interactive = interactive;

[cfg, artifact_EOG] = ft_artifact_zvalue(cfg,data);

% Remove Artifacts from data

clean_data = ft_rejectartifact(cfg,data);

end

% clean_data = all_data; % bypass rejection
% 
% % Concatonate artifact-free trials
% 
% dataConcat = [];
% timeConcat = [];
% 
% for i = 1:size(clean_data.trial,2)
%     dataConcat = [dataConcat, clean_data.trial{i}];
%     if i == 1
%         nextTime = clean_data.time{i};
%     elseif i > 1
%         nextTime = clean_data.time{i} + max(timeConcat);
%     end
%     timeConcat = [timeConcat, nextTime];
% end
% 
% dataConcat(isnan(dataConcat)) = 0; % turn nans into zero
% 
% concatData = [];
% 
% % concatData.hdr = clean_data.hdr;
% concatData.label = clean_data.label;
% concatData.fsample = clean_data.fsample;
% concatData.trial = {dataConcat};
% concatData.time = {timeConcat};
% %     concatData.dimord = 'chan_time';
% %     concatData.sampleinfo = [1 104381];
% %     concatData.elec = ft_read_sens(fileID);
% %     concatData.trialinfo = clean_data.trialinfo;
% %     concatData.sampleinfo = clean_data.sampleinfo;
% %     concatData.cfg = clean_data.cfg;
% 
% % padData = ft_preproc_padding(concatData,'zero',250); % avoid padding error?