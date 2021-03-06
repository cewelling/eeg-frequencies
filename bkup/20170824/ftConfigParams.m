function [ all_data ] = ftConfigParams( parName, runName, EEGfile, numTrials, customDur, date, event, lowBound, highBound)
% lowBound and highBound are the bounds of the band pass filter
%UNTITLED15 Summary of this function goes here
%   Returns pre-processed data with noisy electrodes eliminated (all_data),


% load parameters
analysisParams

if strcmp(parName, 'cumulus25')
    elecsUsed = 33:64;
else
    elecsUsed = 1:32;
end

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

% these eventvalues were defined in runScript, which wrote bytes to EEG comp port
% when participant pressed "up" to begin trial

% Mistakes-----------------------------------------------------------------
% SHOULD NOT NEED THIS (201 and 202 signals should not exist)
% % Did not save first two trials
% if ~isempty(strfind(EEGfile, 'cumulus05_marzRival1'))
%     cfg_trldef.trialdef.eventvalue = 203:206;
% end

%--------------------------------------------------------------------------

try
    cfg_trldef = ft_definetrial(cfg_trldef);
catch define_trial_error
    cfg_trldef.trialdef.eventvalue
    cfg_trldef.dataset
    rethrow(define_trial_error);
end

%% Pre-pre-processing (Identify noisy electrodes)

if ~exist(['pre-processing/badElecs/' parName '_' runName '.mat'], 'file')
    cfg_preproc = cfg_trldef;
    cfg_preproc.channel = elecsUsed;
    cfg_preproc.continuous = 'yes';
    cfg_preproc.demean    = 'no'; % perform baseline correction
    
    % for a few trials, baseline window comes after signal
    if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
            || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
            || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
        %cfg_preproc.baselinewindow = [0 1.5];
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
    
    %% Remove trials that aren't real (in pre-pre-processing file)
    % A mistake in our runScripts was sending some misplaced trial-marking signals to
    % the EEG computer
    
    if strcmp(event, 'trial')
        thisData = checkNoise_data;
        
        % stratus138
        if ~isempty(strfind(EEGfile, 'stratus138_dartRival1')) || ~isempty(strfind(EEGfile, 'stratus138_dartSim3')) ...
                || ~isempty(strfind(EEGfile, 'stratus138_marzRival1'))
            thisData = removeTrials(thisData, numTrials, 3);
        elseif ~isempty(strfind(EEGfile, 'stratus138_dart')) || ~isempty(strfind(EEGfile, 'stratus138_marzRival2'))
            thisData = removeTrials(thisData, numTrials, 2, 4);
        elseif ~isempty(strfind(EEGfile, 'stratus138'))
            thisData = removeTrials(thisData, numTrials, 2);
        end
        
        % cumulus15
        if ~isempty(strfind(EEGfile, 'cumulus15_dartSim1')) || ~isempty(strfind(EEGfile, 'cumulus15_dartSim2')) ...
                || ~isempty(strfind(EEGfile, 'cumulus15_marzRival3'))
            thisData = removeTrials(thisData, numTrials, 3);
        elseif ~isempty(strfind(EEGfile, 'cumulus15_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus15_marzRival'))
            thisData = removeTrials(thisData, numTrials, 2, 4);
        elseif ~isempty(strfind(EEGfile, 'cumulus15'))
            thisData = removeTrials(thisData, numTrials, 2);
        end
        
        % cumulus05
        if ~isempty(strfind(EEGfile, 'cumulus05_dartRival2')) || ~isempty(strfind(EEGfile, 'cumulus05_dartSim2'))
            thisData = removeTrials(thisData, numTrials, 2, 4);
        elseif ~isempty(strfind(EEGfile, 'cumulus05_dartSim3')) || ~isempty(strfind(EEGfile, 'cumulus05_marzRival2')) ...
                || ~isempty(strfind(EEGfile, 'cumulus05_marzRival3'))
            thisData = removeTrials(thisData, numTrials, 2);
        end
        
        % cumulus06
        if ~isempty(strfind(EEGfile, 'cumulus06_marzRival1'))
            thisData = removeTrials(thisData, numTrials, 3);
        elseif ~isempty(strfind(EEGfile, 'cumulus06_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus06_dartSim1')) ...
                || ~isempty(strfind(EEGfile, 'cumulus06_dartSim3')) || ~isempty(strfind(EEGfile, 'cumulus06_marzRival'))
            thisData = removeTrials(thisData, numTrials, 2, 4);
        elseif ~isempty(strfind(EEGfile, 'cumulus06'))
            thisData = removeTrials(thisData, numTrials, 2);
        end
        
        % cumulus28
        if ~isempty(strfind(EEGfile, 'cumulus28_dartRival1'))
            thisData = removeTrials(thisData, numTrials, 1);
        end
        
        checkNoise_data = thisData;
    end
    
    %% Create histograms for noisy electrode detection
    kurMatrix = [];
    meanMatrix = [];
    varMatrix = [];
    for iTrial = 1:length(checkNoise_data.trial)
        trialData = checkNoise_data.trial{iTrial};
        %for elec = 1:size(trialData,1)
        for elec = 1:32
            kur = kurtosis(trialData(elec,:));
            mn = mean(trialData(elec,:));
            vr = var(trialData(elec,:));
            kurMatrix = [kurMatrix kur];
            meanMatrix = [meanMatrix mn];
            varMatrix = [varMatrix vr];
        end
    end
    
    figure
    hist(kurMatrix, 250);
    title([runName ' kurtosis'])
    drawnow
    adjustHist(kurMatrix)
    fprintf('choose min \n');
    [kurMin,~] = ginput(1);
    fprintf('kurtosis min selected \n');
    fprintf('choose max \n');
    [kurMax,~] = ginput(1);
    fprintf('kurtosis max selected \n');
    close gcf
    
    figure
    hist(meanMatrix, 250)
    title([runName ' mean'])
    drawnow
    adjustHist(meanMatrix)
    fprintf('choose min \n');
    [meanMin,~] = ginput(1);
    fprintf('mean min selected \n');
    fprintf('choose max \n');
    [meanMax,~] = ginput(1);
    fprintf('mean max selected \n');
    close gcf
    
    figure
    hist(varMatrix, 250)
    title([runName ' variance'])
    drawnow
    adjustHist(varMatrix)
    fprintf('choose min \n');
    [varMin,~] = ginput(1);
    fprintf('variance min selected \n');
    fprintf('choose max \n');
    [varMax,~] = ginput(1);
    fprintf('variance max selected \n');
    close gcf
    
    %% Find noisy electrodes
    badElecs = {};
    
    for iTrial = 1:length(checkNoise_data.trial)
        trialData = checkNoise_data.trial{iTrial};
        badCount = 1;
        %         for iElec = 1:size(trialData,1)
        for iElec = 1:32
            kur = kurtosis(trialData(iElec,:));
            vr = var(trialData(iElec,:));
            mn = mean(trialData(iElec,:));
            if kur < kurMin || kur > kurMax ...
                    || vr < varMin || vr > varMax ...
                    || mn < meanMin || mn > meanMax
                badElecs{iTrial}(badCount).num = iElec;
                badElecs{iTrial}(badCount).kur = kur;
                badElecs{iTrial}(badCount).var = vr;
                badElecs{iTrial}(badCount).mean = mn;
                badCount = badCount + 1;
            end
        end
    end
    save(['pre-processing/badElecs/' parName '_' runName], 'badElecs');
    
    %     cfg_rejectvis = [];
    %     cfg_rejectvis.method = 'summary';
    %     cfg_rejectvis.keepchannel = 'nan';
    %     pruned_data = ft_rejectvisual(cfg_rejectvis,checkNoise_data);
    %
    %     % store noisy electrodes
    %     badElecs = [];
    %     trialData = pruned_data.trial{1}; % all trials have the same electrodes eliminated, only need to check one
    %     badCount = 1;
    %     for iElec = 1:size(trialData,1);
    %         if isnan(trialData(iElec, 1))
    %             badElecs(badCount).num = iElec;
    %             badCount = badCount+1;
    %         end
    %     end
    %     save(['badElecs/' parName '_' runName], 'badElecs');
    
else
    load(['pre-processing/badElecs/' parName '_' runName '.mat']);
end

%% Preprocessing

cfg_preproc = cfg_trldef;
cfg_preproc.channel = elecsUsed;
cfg_preproc.continuous = 'yes';
cfg_preproc.demean    = 'no'; % perform baseline correction

% for a few trials, baseline window comes after signal
if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
    %cfg_preproc.baselinewindow = [0 1.5];
end

% Started adding a baseline window before signal
date = datetime(date, 'Format','yyyy-mm-dd');
if date > datetime('2016-12-02','Format','yyyy-mm-dd');
    %cfg_preproc.baselinewindow = [-1.5 0];
end

cfg_preproc.detrend = 'yes'; % remove linear trend from the data

% re-reference
cfg_preproc.reref = 'yes'; %%%'yes'; 
refElecs = 1:32;
% for i = 1:length(refElecs)
%     refElecs{i} = num2str(refElecs{i});
% end
% refElecs = [refElecs]; %, 'veog'; 'heog'; 'lm'; 'rm'}; % veog and heog are mastoid channels, lm and rm are eyes

% preclude noisy electrodes from being included as a reference

% get list of all electrodes that were noisy for at least one trial during
% this run
%**************************************************************************
% runBadElecs = [];
% for iTrial = 1:length(badElecs)
%     if isempty(badElecs{iTrial})
%         runBadElecs = union(runBadElecs, []);
%     % if > half of electrodes were bad in this trial, assume globally bad trial; don't include electrodes in bad elec list
%     elseif length(badElecs{iTrial}) < 15
%         runBadElecs = union(runBadElecs, extractfield(badElecs{iTrial},'num'));
%     end
% end
% 
% % Get labels and indices for bad electrodes
% % goodIndices = 1:length(refElecs);
% for i = 1:length(runBadElecs)
%     %     if runBadElecs(i) < 65
%     %         badElec = num2str(runBadElecs(i));
%     %     elseif runBadElecs(i) == 65
%     %         badElec = 'veog';
%     %     elseif runBadElecs(i) == 66
%     %         badElec = 'heog';
%     %     elseif runBadElecs(i) == 67
%     %         badElec = 'lm';
%     %     elseif runBadElecs(i) == 68
%     %         badElec = 'rm';
%     %     else
%     %         badElec = [];
%     %     end
%     %     badIndex = find(strcmp(refElecs, badElec));
%     %     if ~isempty(badIndex)
%     %         goodIndices = goodIndices(goodIndices ~= badIndex);
%     % reference only good electrodes
%     refElecs = refElecs(refElecs ~= runBadElecs(i));
% end
%**************************************************************************
%end

% reference only good electrodes
%refElecs = refElecs(goodIndices);

cfg_preproc.refchannel = refElecs; % common average or median reference
%cfg_preproc.refchannel = 'all'; % common average reference

cfg_preproc.refmethod = 'avg'; %'avg'; %'avg'; % median %cer 11/16
% cfg_preproc.implicitref = 'M1';            % the implicit (non-recorded) reference channel is added to the data representation
% cfg_preproc.refchannel     = {'M1', 'M2'}; % the average of these channels is used as the new reference, note that channel '53' corresponds to the right mastoid (M2)

%cfg_preproc.hpfilter = 'yes'; % highpass filter
%cfg_preproc.hpfreq = 2; % highpass frequency in hertz
cfg_preproc.bpfilter = 'yes'; % bandpass filter?
cfg_preproc.bpfreq = [lowBound highBound]; %[4 10]; %[2 70]; % bandpass frequencies
%cfg_preproc.bpfreq = [1 30]; % bandpass frequencies

cfg_preproc.bsfilter = 'yes'; % bandstop filter?
cfg_preproc.bsfreq = [59 61]; % bandstop frequencies (US = 60)

all_data = ft_preprocessing(cfg_preproc);

%% Remove trials that aren't real
% A mistake in our runScripts was sending some misplaced trial-marking signals to
% the EEG computer

if strcmp(event, 'trial')
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
    
    % cumulus28
    if ~isempty(strfind(EEGfile, 'cumulus28_dartRival1'))
        all_data = removeTrials(all_data, numTrials, 1);
    end
end

%% Remove bad electrodes (on a trial-by-trial basis)
if strcmp(event, 'trial')
    badCount = zeros(32,1);
    for iTrial = 1:length(badElecs)
        for iElec = 1:length(badElecs{iTrial})
            badElec = badElecs{iTrial}(iElec).num;
            % replace bad electrode with nans
            all_data.trial{iTrial}(badElec,:) = nan(1,size(all_data.trial{iTrial}(badElec,:),2));
            
            badCount(badElec) = badCount(badElec) + 1;
        end
    end
    
    % if electrode was noisy in more than three trials, remove it altogether
    for i = 1:length(badCount)
        if badCount(i) > 3
            for iTrial = 1:length(badElecs)
                all_data.trial{iTrial}(i,:) = nan(1,size(all_data.trial{iTrial}(i,:),2));
            end
        end
    end
end

%% Test electrodes (for finding a metric for noisy electrodes)
elecsChosen = []; % include elecs to be analyzed, plotted
testElecs = {};
for iTrial = 1:length(all_data.trial)
    trialData = all_data.trial{iTrial};
    count = 1;
    for elec = elecsChosen;
        kur = kurtosis(trialData(elec,:));
        vr = var(trialData(elec,:));
        mn = mean(trialData(elec,:));
        testElecs{iTrial}(count).num = elec;
        testElecs{iTrial}(count).kur = kur;
        testElecs{iTrial}(count).var = vr;
        testElecs{iTrial}(count).mean = mn;
        count = count + 1;
    end
end
save(['pre-processing/testElecs/' parName '_' runName], 'testElecs');

hhh
end

