function [ rls_data, rls_time ] = runRLS(parName, runName, date, EEGfile)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
analysisParams

%% Trial definition and pre-processing

all_data = ftConfigParams(parName, runName, EEGfile, numTrials, trialDur, date);

 %% RLS Analysis

 rls_data = [];
 msglength = 0;
 analElecs = load(['highSNRelecs/' parName '_' runName '_'  num2str(numElecs) 'elecs.mat']);
 analElecs = analElecs.maxIndices;
 for iTrial = 1:size(all_data.trial, 2); % # trials recorded in EEG (should = numTrials)
     
     fprintf(repmat('\b',1,msglength));
     msg = sprintf('Running RLS - Processing Trial No %d\n', iTrial);
     fprintf(msg);
     msglength = numel(msg);
     
     for iFreq = 1:length(stimFreqs) % size([stimFreqs imFreqs])
         currFreq = stimFreqs(iFreq);
         window = 1000/currFreq*5;
         window = timeWindow;
         sampleRate = 512;
         
         % Use the RLS filter
         %[Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(electrodeSet.nums,:))',currFreq,[],1,sampleRate);
         [Hcn1,Hsn1] = AdaptiveRLSfast(window,(all_data.trial{iTrial}(analElecs,:))',currFreq,[],1,sampleRate);
         
         cplx = complex(Hcn1,Hsn1);
         envelope = abs(cplx);
         envelope = envelope';
         
         for iElec = 1:length(electrodeSet.nums)
             % smooth RLS amplitudes
             envelope(iElec, :) = smoothts(envelope(iElec,:), 'g', smoothingWin,smoothingStd);
             %envelope(iElec, :) = smooth(envelope(iElec,:), smoothingPts);
             
             % eliminate first window of data, shift all remaining data over
             pad = round((window/1000)*sampleRate);
             shiftedEnv(iElec, :) = [envelope(iElec, (pad+1):end) (mean(envelope((pad+1):end))*ones(1,pad))];
             % pads end with average of all the data
         end
         
         % average across electrodes
         shiftedEnv = mean(shiftedEnv, 1);
            
         % store time vs. smoothed RLS amplitude
         thisTrial = iTrial;
         if ~isempty(strfind(EEGfile, 'cumulus05_marzRival1'))
             % account for loss of first two trials of this run
             thisTrial = thisTrial + 2;
         end
         rls_data(iFreq).amp {thisTrial} = shiftedEnv;
         rls_data(iFreq).electrodes {thisTrial} = electrodeSet.name;
     end
 end
 
 rls_time = all_data.time{1}; % time vector should be the same for all trials
 
 
 %% Plot the time series, overlay perceptual epochs
 
 if strcmp(RLSplotOrNot, 'yes')
     for iTrial = 1:size(all_data.trial, 2)
         
         trialTime = rls_time; 
         
         % Average across electrodes
         lowBand = mean(rls_data(1).amp{iTrial}, 1);
         lowBandError = ste(rls_data(1).amp{iTrial});
         
         highBand = mean(rls_data(2).amp{iTrial}, 1);
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

