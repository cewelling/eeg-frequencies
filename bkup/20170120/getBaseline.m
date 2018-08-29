function [baselines] = getBaseline( parName, runName)
% Find baseline by averaging periods without stimuli on screen during a 
% simulation run

% load parameters
analysisParams

% Get EEG file info
EEGfiles = dir([eegDir parName '*']);
% ...if the participant exists...
if isempty(EEGfiles)
    error('Participant does not exist - cannot get baseline');
end
date = strtok(EEGfiles(1).date);
date = datestr(date, dateFormat);
EEGfile = [eegDir parName '_' runName '_' date '.bdf'];

%% baseline window is in different places for different participants / runs
dateT = datetime(date, 'Format','yyyy-MM-dd');
% baseline window is after signal that starts trial
% if simBase == 1
%     type = 'simBase';
if ~isempty(strfind(EEGfile, 'cumulus10_dartRival1')) || ~isempty(strfind(EEGfile, 'cumulus10_dartRival2')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartRival3')) || ~isempty(strfind(EEGfile, 'cumulus10_dartSim1')) ...
        || ~isempty(strfind(EEGfile, 'cumulus10_dartSim2')) || ~isempty(strfind(EEGfile, 'cumulus13'))
    type = 'postStim'; %'preStim2';
    % Started adding a baseline window before signal
elseif dateT > datetime('2016-12-02','Format','yyyy-MM-dd');
    type = 'preStim'; %'preStim';
else
    type = 'trial';%'postStim';
end

% load key press data file
keyPressData = load([keyPressDir parName '_' runName '_' date '.txt']);


%% Pre-processing
all_data = ftConfigParams(parName, runName, EEGfile, numTrials, 1.5, date, type);

%% RLS Analysis

 rls_data = [];
 msglength = 0;
 
 % get electrodes for this run
 if strfind(runName, 'dart')
     analElecs = load(['highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
 else
     analElecs = load(['highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
 end
 
 if isempty(analElecs.elecs) % SNR too low in this run
     error('SNR too low in this run for use as baseline');
 end
 
 analElecs = analElecs.elecs(1,:);
 
 % Run RLS analysis to get baseline amplitudes
 for iTrial = 1:size(all_data.trial, 2); % # trials recorded in EEG (should = numTrials)
     
     fprintf(repmat('\b',1,msglength));
     msg = sprintf('Running RLS - Processing Trial No %d\n', iTrial);
     fprintf(msg);
     msglength = numel(msg);
     
     %for iFreq = 1:length(stimFreqs) % size([stimFreqs imFreqs])
     for stimFreq = 1:length(stimFreqs) % size([stimFreqs imFreqs])
         % determine frequencies for baseline calculation
         rHalf = 1.5; % range: stimFreq +/- freqRange
         numFreqs = 10;
         step = rHalf / numFreqs * 2;
         bLineFreqs = [stimFreq-rHalf:step:stimFreq - step stimFreq+step:step:stimFreq+rHalf];
         
         for iFreq = 1:length(bLineFreqs) % size([stimFreqs imFreqs])
             %currFreq = stimFreqs(iFreq);
             currFreq = bLineFreqs(iFreq);
             
             window =  sampRate/currFreq*10;
             
             % Use the RLS filter
             %[Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(electrodeSet.nums,:))',currFreq,[],1,sampleRate);
             [Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(analElecs,:))',currFreq,[],1,sampRate);
             
             cplx = complex(Hcn1,Hsn1);
             envelope = abs(cplx);
             envelope = envelope';
             
             % average across electrodes
             envelope = nanmean(envelope, 1);
             
             % store time vs. smoothed RLS amplitude
             thisTrial = iTrial;
             
             % Account for mistakes
             if ~isempty(strfind(EEGfile, 'cumulus05_marzRival1'))
                 % account for loss of first two trials of this run
                 thisTrial = thisTrial + 2;
             elseif ~isempty(strfind(EEGfile, 'stratus135_dartSim1'))
                 % account for loss of first trial of this run
                 thisTrial = thisTrial + 1;
             elseif ~isempty(strfind(EEGfile, 'cumulus01_dartSim3'))
                 % account for loss of first trial of this run
                 thisTrial = thisTrial +1;
             elseif ~isempty(strfind(EEGfile, 'cumulus01_marzRival3'))
                 % account for loss of first two trials of this run
                 thisTrial = thisTrial +2;
             end
             figure
             plot(envelope)
             title(num2str(bLineFreqs(iFreq)))
             ylim([0 6])
             % store average of RLS amplitudes
             rls_data(iFreq,thisTrial) = mean(envelope(300:end));
         end
     end
 end
 
 baselines = rls_data;
 save(['baselineTest/' parName '_' runName '_' type], 'baselines')
 
 % average across trials
 baselines = mean(rls_data, 2);
end
