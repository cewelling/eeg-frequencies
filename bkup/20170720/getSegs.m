function [ lowSegList, highSegList, labelVect ] = getSegs( parName, date, trainType, domSegProp, mixSegProp, segPoints, channel, freqType, numStates, paramsFile )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
run(paramsFile);

ZS = 0; % z-score?
unsmoothed = 0; % use unsmoothed RLS data?

% space to store lists of segments and labels
lowSegList = [];
highSegList = [];
labelVect = {}; % category vector

%% Iterate through runs to find epoch segments
for iRunType = trainType % just input 1 for dartRival only and 2 for dartSim only
    cRunType = runTypes{iRunType};
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        % Get the appropriate EEG file for this run
        EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
        
        % skip over participants / runs that don't exist
        if ~exist(EEGfile, 'file')
            continue;
        end
        
        %% Get RLS data
%         if unsmoothed % unsmoothed
%             if ~strcmp(freqType, 'principle')
%                 error('Script not yet set up for non-principle frequencies with unsmoothed RLS data');
%             end
%             [rls_data, rls_time] = runRLS_noSmoothing(parName, runName, date, EEGfile);
%         else % smoothed
%             if strcmp(freqType, 'principle') % principle stimulation frequencies
%                 if exist(['rls_data/smoothedRLS/' parName '_' runName '.mat'], 'file')
%                     load(['rls_data/smoothedRLS/' parName '_' runName '.mat']);
%                 else
%                     if analFreqs == stimFreqs
%                         [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
%                     else
%                         error('In analysisParams.m, set analFreqs = stimFreqs, then re-run');
%                     end
%                 end
%             elseif strcmp(freqType, 'harmonics') % harmonics of stimulation frequencies
%                 if exist(['rls_data/smoothedRLS/' parName '_' runName '_harmonics.mat'], 'file')
%                     load(['rls_data/smoothedRLS/' parName '_' runName '_harmonics.mat']);
%                 else
%                     if analFreqs == harFreqs
%                         [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
%                     else
%                         error('In analysisParams.m, set analFreqs = harFreqs, then re-run');
%                     end
%                 end
%             elseif strcmp(freqType, 'im') % intermodulation frequencies
%                 if exist(['rls_data/smoothedRLS/' parName '_' runName '_im.mat'], 'file')
%                     load(['rls_data/smoothedRLS/' parName '_' runName '_im.mat']);
%                 else
%                     if analFreqs == imFreqs
%                         [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
%                     else
%                         error('In analysisParams.m, set analFreqs = imFreqs, then re-run');
%                     end
%                 end
%             end
%         end
        
        %% get list of participant's perceptual epochs
        if iRunType == 2 % use sim trials
            if ~exist(['simSchedules/' parName '_' runName '.mat'], 'file')
                [~, simSchedule] = getEpochs(parName, runName, date);
            else
                load(['simSchedules/' parName '_' runName '.mat'])
            end
            keyPressRef = simSchedule;
        elseif iRunType == 1 % use rivalry trials
            if ~exist(['epochs/' parName '_' runName '.mat'], 'file')
                [epochs, ~] = getEpochs(parName, runName, date);
            else
                load(['epochs/' parName '_' runName '.mat'])
            end
            keyPressRef = epochs;
        end
        
        % Ensure that individual channel rls analysis has been completed for this run
        thisRunFiles = dir(['rls_data/perChannel/' parName '_' runName '*' num2str(channel) '_' freqType '.mat']);
        if isempty(thisRunFiles)
            error(['Run RLS analysis from analysisController with electrode ' num2str(channel) ' and ' freqType ' selected']);   
        end
        
        %% Handle each trial in the RLS data
        for iTrial = 1:numTrials
            
            % Load RLS data
            if strcmp(freqType, 'principle')
                if exist(['rls_data/perChannel/' parName '_' runName '_trial' num2str(iTrial) '_freq1_' num2str(channel) '.mat'], 'file')
                    load(['rls_data/perChannel/' parName '_' runName '_trial' num2str(iTrial) '_freq1_' num2str(channel) '.mat']);
                    lowBand = thisElecRLS;
                    load(['rls_data/perChannel/' parName '_' runName '_trial' num2str(iTrial) '_freq2_' num2str(channel) '.mat']);
                    highBand = thisElecRLS;
                else
                    continue; % account for trials that weren't recorded properly (mistakes)
                end
            else
                if exist(['rls_data/perChannel/' parName '_' runName '_trial' num2str(iTrial) '_freq1_' num2str(channel) '_' freqType '.mat'], 'file')
                    load(['rls_data/perChannel/' parName '_' runName '_trial' num2str(iTrial) '_freq1_' num2str(channel) '_' freqType '.mat']);
                    lowBand = thisElecRLS;
                    load(['rls_data/perChannel/' parName '_' runName '_trial' num2str(iTrial) '_freq2_' num2str(channel) '_' freqType '.mat']);
                    highBand = thisElecRLS;
                else
                    continue; % account for trials that weren't recorded properly (mistakes)
                end
            end
            
            if ZS
                % Z-score
                lowBand(~isnan(lowBand)) = zscore(lowBand(~isnan(lowBand))); % ignore NaNs when z-scoring
                highBand(~isnan(highBand)) = zscore(highBand(~isnan(highBand)));
            else
                % Just demean
                lowBand = lowBand - nanmean(lowBand);
                highBand = highBand - nanmean(highBand);
            end
            
            % is the SNR of this trial high enough (in both bands)?
            load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' electrodeSet.name '.mat']);
            %if nanmean(lFreqSNRs(2,lFreqSNRs(1,:) == channel)) < minSNR || ...
            %nanmean(hFreqSNRs(2,hFreqSNRs(1,:) == channel)) < minSNR
            if nanmean(lFreqSNRs(2,lFreqSNRs(1,:) == 29)) < minSNR || ...
                    nanmean(hFreqSNRs(2,hFreqSNRs(1,:) == 29)) < minSNR
                continue;
            end
            
            % Get perceptual epochs for this trial
            trialEpochs = keyPressRef(keyPressRef(:,1) == iTrial,:);
            
            %% Plot RLS traces to visualize segments
            if strcmp(segPlotOrNot, 'yes')
                
                % load RLS time (pick a random smoothed RLS file)
                load('rls_data/smoothedRLS/cumulus01_dartRival1_principles.mat')
                
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
            
            %% partition RLS data into epochs based on sim schedule
            for iEpoch = 1:size(trialEpochs,1)
                
                % eliminate first mixed epoch
                if iEpoch == 1 && trialEpochs(iEpoch, 5) == 0
                    continue;
                end
                
                eStart = trialEpochs(iEpoch,2); % in seconds
                eEnd = trialEpochs(iEpoch,3); % in seconds
                
                % minimum epoch length
                if eEnd - eStart < minEpochDur % minEpochDur = 0.5 (in sec)
                    continue;
                end
                
                % Position of segment...
                
                %eMid = eStart + segMiddle; %eStart + ((eEnd - eStart)/ 2); %eEnd - segMiddle; %**eStart + segMiddle**; %eEnd; %eStart + ((eEnd - eStart)/ 2); % in seconds
                %eMid = eEnd - segMiddle;
                %eMid = eStart + ((eEnd - eStart)/ 2);
                % ... for mixed epoch
                if trialEpochs(iEpoch, 5) == 0
                    eMid = eStart + ((eEnd - eStart) * mixSegProp);
                    if eEnd - eMid < segPoints / sampRate / 2
                        eMid = eEnd - segPoints / sampRate / 2;
                    elseif eMid - eStart < segPoints / sampRate / 2
                        eMid = eStart + segPoints / sampRate / 2;
                    end
                % ... for dominant epoch
                else
                    %eMid = eEnd - domSegMid; % UNCOMMENT AFTER TESTING
                    %eMid = eStart + ((eEnd - eStart) * 0.5); % TEMP FOR TESTING
                    eMid = eStart + ((eEnd - eStart) * domSegProp);
                    if eEnd - eMid < segPoints / sampRate / 2
                        eMid = eEnd - segPoints / sampRate / 2;
                    elseif eMid - eStart < segPoints / sampRate / 2
                        eMid = eStart + segPoints / sampRate / 2;
                    end
                end
%                 if eEnd - eStart < segMiddle
%                     eMid = eEnd;
%                 end
                startI = round(eStart*sampRate);
                eMid = eMid - 1; % account for cutting 1 second off beginning of RLS data
                midI = round(eMid*sampRate); % in data points
                
                % NO LATENCY WHEN USING SIM SCHEDULE
                % account for latency
                % midI_lat = midI - latency;
                
                segStart = midI - round(segPoints/2);
                %segStart = startI;
                segEnd = segStart + segPoints - 1;
                
                
                % Segment includes points within first 1 or last 1.5 seconds of trial; eliminate it
                if segStart <= sampRate
                    continue;
                elseif segEnd >= sampRate*(trialDur-1.5)
                    continue;
                end
                
                % get segment
                seg = segStart:segEnd;
                lowSeg = lowBand(seg);
                highSeg = highBand(seg);
                
                % Bin the data according to the given bin interval
                %                 leftOver = mod(segPoints, binPoints);
                %                 for currPoint = ceil(leftOver/2)+binPoints:binPoints:segPoints - floor(leftOver/2)
                %                     lowBinVal = nanmean(lowSeg(currPoint - binPoints:currPoint - 1));
                %                     highBinVal = nanmean(highSeg(currPoint - binPoints:currPoint - 1));
                %                     lowFreqVect = [lowFreqVect; lowBinVal];
                %                     highFreqVect = [highFreqVect; highBinVal];
                %                 end
                
                % collect segments
                lowSegList = [lowSegList; lowSeg];
                highSegList = [highSegList; highSeg];
                
                % Shade segment in RLS plot and store its label
                if strcmp(segPlotOrNot, 'yes')
                    %horiz = [seg(1)/sampRate seg(end)/sampRate];
                    horiz = [eStart - 1 eEnd - 1]; % account for cutting 1 second off beginning of EEG data
                    vertUp = [uShadeLim uShadeLim];
                    vertDown = [lShadeLim lShadeLim];
                end
                
                if trialEpochs(iEpoch,5) > 0.5
                    if strcmp(segPlotOrNot, 'yes')
                        area(horiz,vertUp,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
                        area(horiz,vertDown,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
                    end
                    if numStates == 3
                        labelVect = [labelVect; 'hi'];
                    elseif numStates == 2
                        labelVect = [labelVect; 'dom'];
                    end
                    
                elseif trialEpochs(iEpoch,5) < -0.5
                    if strcmp(segPlotOrNot, 'yes')
                        area(horiz,vertUp,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
                        area(horiz,vertDown,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
                    end
                    if numStates == 3
                        labelVect = [labelVect; 'lo'];
                    elseif numStates == 2
                        labelVect = [labelVect; 'dom'];
                    end
                else
                    if strcmp(segPlotOrNot, 'yes')
                        area(horiz,vertUp,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
                        area(horiz,vertDown,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
                    end
                    labelVect = [labelVect; 'mi'];
                end
            end
        end
    end
end

end

