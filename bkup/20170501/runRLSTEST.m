function [ rls_data, rls_time ] = runRLSTEST(parName, runName, date, EEGfile)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
analysisParams

 %% Get stored analysis electrodes 

 rls_data = [];
 msglength = 0;
 
 % Get stored analysis electrodes 
 if strfind(runName, 'dart')
     analElecs = load(['highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
 else
     analElecs = load(['highSNRelecs/' parName '_marz_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
 end
 
 if isempty(analElecs.elecs) % SNR too low in this run (not an exhaustive snr gate)
     rls_data = [];
     rls_time = [];
     return;
 end
 
 analElecs = analElecs.elecs(1,:);
 
 for iFreq = 1:length(analFreqs) % size([stimFreqs imFreqs])
     
     %% Trial definition and pre-processing
     
     if exist(['preproc/' parName '_' runName '_freq' num2str(iFreq) '.mat'], 'file')
         load(['preproc/' parName '_' runName '_freq' num2str(iFreq) '.mat'])
     else
         % Ensure that path is in correct order
         path(genpath('/Users/acspiegel/Documents/MATLAB/Toolboxes/FieldTrip'), path)
         
         all_data = ftConfigParams(parName, runName, EEGfile, numTrials, trialDur, date, [], analFreqs(iFreq) - 2, analFreqs(iFreq) + 2);
         save(['preproc/' parName '_' runName '_freq' num2str(iFreq)], 'all_data');
     end
     
     for iTrial = 1:size(all_data.trial, 2); % # trials recorded in EEG (should = numTrials)
         if exist(['unsmoothedRLS/' parName '_' runName '_trial' num2str(iTrial) '_freq' num2str(iFreq) '.mat'], 'file')
             load(['unsmoothedRLS/' parName '_' runName '_trial' num2str(iTrial) '_freq' num2str(iFreq) '.mat'])
         else
             fprintf(repmat('\b',1,msglength));
             msg = sprintf('Running RLS - Processing Trial No %d\n', iTrial);
             fprintf(msg);
             msglength = numel(msg);
             
             %for iFreq = 1:length(stimFreqs) % size([stimFreqs imFreqs])
             currFreq = analFreqs(iFreq);
             window = sampRate/currFreq*10; %1000/currFreq*5; (in timePoints)
             
             % Use the RLS filter
             %[Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(electrodeSet.nums,:))',currFreq,[],1,sampleRate);
             [Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(electrodeSets.nums(analElecs),:))',currFreq,[],1,sampRate);
             
             cplx = complex(Hcn1,Hsn1);
             envelope = abs(cplx);
             envelope = envelope';
             
             save(['unsmoothedRLS/' parName '_' runName '_trial' num2str(iTrial) '_freq' num2str(iFreq)], 'envelope');
         end
         
         for iElec = 1:length(analElecs)
             % smooth RLS amplitudes
             %envelope(iElec, :) = HRsmoothing(envelope(iElec, :), 'gaussian', 151, smooth_sd, 0);
             %envelope(iElec, :) = smoothts(envelope(iElec,:), 'g', smoothingWin,smoothingStd);
             %envelope(iElec, :) = smooth(envelope(iElec,:), smoothingPts);
             
             % eliminate first second of data, shift all remaining data over
             pad = round((cutoffWindow/1000)*sampRate);
             shiftedEnv(iElec, :) = [envelope(iElec, (pad+1):end) (nanmean(envelope((pad+1):end))*ones(1,pad))];
             % pads end with average of all the data
         end
         
         % average across electrodes
         shiftedEnv = nanmean(shiftedEnv, 1);
         
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
         rls_data(iFreq).amp {thisTrial} = shiftedEnv;
         rls_data(iFreq).electrodes {thisTrial} = electrodeSet.name;
     end
 end
 
 rls_time = all_data.time{1}; % time vector should be the same for all trials
 
 % save RLS data
 %  if analFreqs == stimFreqs
%      save(['rls_data/' parName '_' runName], 'rls_data', 'rls_time');
%  elseif analFreqs == harFreqs
%      save(['rls_data/' parName '_' runName '_harmonics'], 'rls_data', 'rls_time');
%  end
 
 %% Plot the time series, overlay perceptual epochs
 
 if strcmp(RLSplotOrNot, 'yes')
     for iTrial = 1:size(all_data.trial, 2)
         
         trialTime = rls_time;
         
         % Average across electrodes
         lowBand = nanmean(rls_data(1).amp{iTrial}, 1);
         lowBandError = ste(rls_data(1).amp{iTrial});
         
         highBand = nanmean(rls_data(2).amp{iTrial}, 1);
         highBandError = ste(rls_data(2).amp{iTrial});
         
         % get list of participant's perceptual epochs
         [epochs, ~] = getEpochs(parName, runName, date);
         
         %plotting
         trialEpochs = epochs(epochs(:,1) == iTrial,:);
         
         figure
         title(['RLS: ', parName ',' runName ', trial ' num2str(iTrial)])
         hold on
         %         mseb(trialTime,[lowBand; highBand],[lowBandError; highBandError],[],1)
         plot(trialTime,lowBand,'b','linewidth',2)
         plot(trialTime,highBand,'r','linewidth',2)
         
         % find highest value that will be plotted
         checkMax = [max(lowBand+lowBandError) max(highBand+highBandError)];
         shadeLim = ceil(max(checkMax)); % round up
         
         for iEpoch = 1:size(trialEpochs,1)
             
             horiz = [trialEpochs(iEpoch,2) trialEpochs(iEpoch,3)];
             vert = [shadeLim shadeLim];
             
             if trialEpochs(iEpoch,5) > 0.5
                 area(horiz,vert,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
             elseif trialEpochs(iEpoch,5) < -0.5
                 area(horiz,vert,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
             else
                 area(horiz,vert,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
             end
         end
         
     end
 end
 
end

