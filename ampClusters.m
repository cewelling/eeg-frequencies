function ampClusters( parName, paramsFlag )
% ampClusters(parName) collects amplitudes of demeaned RLS traces (at every
% time point in segments of participant-reported dominant, mixed, and
% suppressed percepts. Saves those amplitudes in the ampClusters folder.
% Plots RLS traces with segments of data used to collect amplitudes
% overlaid. Also visualizes clusters with bar graphs.
%
% This analysis seeks to ask: Do RLS amplitudes cluster according to
% perceptual report? (Dominant vs. Mixed vs. Suppressed?) Saved amplitudes
% are used for group analysis in groupAmpClusters.m.
%
% *** SHIFTS transitions 1 second back to account for 1 second shift of RLS 
% data!!! ***
%
% Called from: analysisController.m
% Dependencies: analysisParams.m

%% Set-up

% Load parameters
analysisParams

% Get EEG file info
EEGfiles = dir([eegDir parName '*']);
% ...if the participant exists...
if isempty(EEGfiles)
    return;
end
date = strtok(EEGfiles(1).date);
date = datestr(date, dateFormat);

%% Iterate through each of the participant's runs and trials

for iRunType = 1:2 %1:length(runTypes) % just input 1 for dartRival only
    cRunType = runTypes{iRunType};
    
    % space to store amplitudes corresponding to each percept
    hMixedAmps = [];
    hCorrAmps = []; % corresponding percept (high frequency if its the high frequency trace)
    hOppAmps = []; % opposite percept (low frequency of its the high frequency trace)
    lMixedAmps = [];
    lCorrAmps = [];
    lOppAmps = [];
    
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        % Get the appropriate EEG file for this run
        EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
        
        % skip over participants / runs that don't exist
        if ~exist(EEGfile, 'file')
            continue;
        end
        
        % Load RLS data
        if smoothing
            if exist(['rls_data/smoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name '.mat'], 'file')
                load(['rls_data/smoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name '.mat']);
            else
                [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile, paramsFlag);
            end
        else
            if exist(['rls_data/unsmoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name '.mat'], 'file')
                load(['rls_data/unsmoothedRLS/' parName '_' runName '_' analFreqLabel '_' electrodeSet.name '.mat']);
            else
                [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile, paramsFlag);
            end
        end
        
        % get list of participant's perceptual epochs
        [epochs, ~] = getEpochs(parName, runName, date, paramsFlag);
        
        % Handle each trial in the RLS data
        for iTrial = 1:size(rls_data(1).amp, 2)
            
            lowBand = (rls_data(1).amp{iTrial});
            highBand = (rls_data(2).amp{iTrial});
            
            % account for trials that weren't recorded properly (mistakes)
            if isempty(lowBand)
                continue;
            end
            
            % is the SNR of this trial high enough?
            load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' electrodeSet.name '.mat']);
            if nanmean(lFreqSNRs(2,lFreqSNRs(1,:)==snrElecs)) < minSNR || ...
                    nanmean(hFreqSNRs(2,hFreqSNRs(1,:)==snrElecs)) < minSNR
                continue;
            end
            
            % Demean the RLS data
            lowBand = lowBand - nanmean(lowBand);
            highBand = highBand - nanmean(highBand);
            
            % Get perceptual epochs for this trial
            trialEpochs = epochs(epochs(:,1) == iTrial,:);
            
            %% Plot RLS traces to visualize segments
            if strcmp(clustSegPlotOrNot, 'yes')
                
                figure
                title(['RLS: ', parName ',' runName ', trial ' num2str(iTrial)])
                hold on
                %         mseb(trialTime,[lowBand; highBand],[lowBandError; highBandError],[],1)
                plot(rls_time,lowBand,'b','linewidth',2)
                plot(rls_time,highBand,'r','linewidth',2)
                
                % find highest value that will be plotted
                checkMax = [max(lowBand) max(highBand)];
                checkMin = [min(lowBand) min(highBand)];
                uShadeLim = ceil(max(checkMax)); % round up
                lShadeLim = floor(min(checkMin)); % round down
            end
            
            %% Collect and save amplitudes from dominant, suppressed, and mixed epochs
            
            for iEpoch = 1:size(trialEpochs,1)
                
                % eliminate first mixed epoch
                if iEpoch == 1 && trialEpochs(iEpoch, 5) == 0
                    continue;
                end
                
                % get start and end of epoch
                % subtract 1 second to account for second cutoff of RLS trace due to stimulus onset
                eStart = trialEpochs(iEpoch,2) - 1; % in seconds 
                eEnd = trialEpochs(iEpoch,3) - 1; % in seconds
                
                % minimum epoch length
                if eEnd - eStart < minEpochDur
                    continue;
                end
                
                startI = round(eStart*sampRate); % in data points
                endI = round(eEnd*sampRate); % in data points
                ePoints = endI - startI; % length of epoch in data points
                
                % account for latency
                startI_lat = startI - latency;
                
                segStart = startI_lat + round(ePoints*(1-segFrac)/2);
                segEnd = segStart + round(ePoints*segFrac);
                
                if segEnd >= sampRate*(trialDur - 1.5)
                    segEnd = sampRate*(trialDur - 1.5) - 1; % eliminate last second of trial
                elseif segStart <= sampRate
                    segStart = sampRate + 1; % eliminate first second of trial
                end
                
                % Epoch is within first 1 or last 1.5 seconds of trial; eliminate it
                if segStart > segEnd
                    continue;
                end
                
                % Shade segments used for analysis in RLS plot and collect amplitudes
                seg = segStart:segEnd;
                
                if strcmp(clustSegPlotOrNot, 'yes')
                    horiz = [seg(1)/sampRate seg(end)/sampRate];
                    vertUp = [uShadeLim uShadeLim];
                    vertDown = [lShadeLim lShadeLim];
                end
                
                if trialEpochs(iEpoch,5) > 0.5
                    if strcmp(clustSegPlotOrNot, 'yes')
                        area(horiz,vertUp,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
                        area(horiz,vertDown,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
                    end
                    lOppAmps = [lOppAmps lowBand(seg)];
                    hCorrAmps = [hCorrAmps highBand(seg)];
                elseif trialEpochs(iEpoch,5) < -0.5
                    if strcmp(clustSegPlotOrNot, 'yes')
                        area(horiz,vertUp,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
                        area(horiz,vertDown,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
                    end
                    lCorrAmps = [lCorrAmps lowBand(seg)];
                    hOppAmps = [hOppAmps highBand(seg)];
                else
                    if strcmp(clustSegPlotOrNot, 'yes')
                        area(horiz,vertUp,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
                        area(horiz,vertDown,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
                    end
                    lMixedAmps = [lMixedAmps lowBand(seg)];
                    hMixedAmps = [hMixedAmps highBand(seg)];
                end
                
            end
        end
    end
    
    save(['ampClusters/' parName '_' cRunType], 'lOppAmps', 'lCorrAmps', 'lMixedAmps', 'hOppAmps', 'hCorrAmps', 'hMixedAmps')
    fprintf('Saving amplitude clusters for group analysis\n')
    
    
    %% Bar plots of amplitudes in each cluster
    
    if strcmp(clustPlotOrNot, 'yes')
        
        % Plot amplitude clusters for low frequency band
        figure
        bar([1 2 3], [nanmean(lCorrAmps) nanmean(lMixedAmps) nanmean(lOppAmps)]);
        hold on
        scatterX = [ones(1, length(lCorrAmps)) (ones(1, length(lMixedAmps)) + 1) (ones(1, length(lOppAmps)) + 2)];
        scatter(scatterX, [lCorrAmps lMixedAmps lOppAmps]);
        title([parName ' ' cRunType ' low frequency band     .'])
        xlabel('     Dominant                                       Mixed                                     Suppressed   ')
        
        % Plot amplitude clusters for high frequency band
        figure
        bar([1 2 3], [nanmean(hCorrAmps) nanmean(hMixedAmps) nanmean(hOppAmps)]);
        hold on
        scatterX = [ones(1, length(hCorrAmps)) (ones(1, length(hMixedAmps)) + 1) (ones(1, length(hOppAmps)) + 2)];
        scatter(scatterX, [hCorrAmps hMixedAmps hOppAmps]);
        title([parName ' ' cRunType ' high frequency band     .'])
        xlabel('     Dominant                                       Mixed                                     Suppressed   ')
        
        histEdges = -3:.2:3;
        figure
        subplot(3,1,1)
        histogram(lCorrAmps, histEdges)
        title('Reported dominant', 'FontWeight', 'normal')
        xlabel('amplitudes', 'FontSize', 16)
        ylabel('# points', 'FontSize', 16)
        subplot(3,1,2)
        histogram(lMixedAmps, histEdges)
        title('Reported mixed', 'FontWeight', 'normal')
        xlabel('amplitudes', 'FontSize', 16)
        ylabel('# points', 'FontSize', 16)
        subplot(3,1,3)
        histogram(lOppAmps, histEdges)
        title('Reported suppressed', 'FontWeight', 'normal')
        xlabel('amplitudes', 'FontSize', 16)
        ylabel('# points', 'FontSize', 16)
        axes('Position',[0 0 1 1],'Visible','off');
        text(0.15,0.98,[parName ' ' cRunType ' low frequency band     .'], 'FontWeight', 'bold')
        
        figure
        subplot(3,1,1)
        histogram(hCorrAmps, histEdges)
        title('Reported dominant', 'FontWeight', 'normal')
        xlabel('amplitudes', 'FontSize', 16)
        ylabel('# points', 'FontSize', 16)
        subplot(3,1,2)
        histogram(hMixedAmps, histEdges)
        title('Reported mixed', 'FontWeight', 'normal')
        xlabel('amplitudes', 'FontSize', 16)
        ylabel('# points', 'FontSize', 16)
        subplot(3,1,3)
        histogram(hOppAmps, histEdges)
        title('Reported suppressed', 'FontWeight', 'normal')
        xlabel('amplitudes', 'FontSize', 16)
        ylabel('# points', 'FontSize', 16)
        axes('Position',[0 0 1 1],'Visible','off');
        text(0.15,0.98,[parName ' ' cRunType ' high frequency band     .'], 'FontWeight', 'bold')
    end
end
end
