%clearvars

% load parameters
%analysisParams

% Get EEG file info
EEGfiles = dir([eegDir parName '*']);
% ...if the participant exists...
if isempty(EEGfiles)
    return;
end
date = strtok(EEGfiles(1).date);
date = datestr(date, dateFormat);

% load or (if not yet set) set common electrodes for each type of run based on SNR
if exist(['pre-processing/highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs.mat'], 'file')
    load(['pre-processing/highSNRelecs/' parName '_darts_' num2str(commonElecs) 'of' num2str(numElecs) 'elecs']);
else
    elecs = setElecs(parName);
end

% Raster data (rows are trials, columns are time points)
highFreqRast = [];
lowFreqRast = []; 

labelVect = {}; % category vector

%% training data set (sim trials)
for iRunType = 2 % just input 2 for dartSim only
    cRunType = runTypes{iRunType};
    for runIndex = runIndices
        runName = [cRunType num2str(runIndex)];
        
        % Get the appropriate EEG file for this run
        EEGfile = [eegDir parName '_' runName '_' date '.bdf'];
        
        % skip over participants / runs that don't exist
        if ~exist(EEGfile, 'file')
            continue;
        end
        
        % Load smoothed RLS data
        if exist(['rls_data/smoothedRLS/' parName '_' runName '.mat'], 'file')
            load(['rls_data/smoothedRLS/' parName '_' runName '.mat']);
        else
            [rls_data, rls_time] = runRLS(parName, runName, date, EEGfile);
        end
        
        % get list of participant's perceptual epochs
        [~, simSchedule] = getEpochs(parName, runName, date);
        
        % Handle each trial in the RLS data
        for iTrial = 1:size(rls_data(1).amp, 2)
            
            % is the SNR of this trial high enough?
            load(['pre-processing/highSNRelecs/' parName '_' runName '_trial' num2str(iTrial) '_' num2str(numElecs) 'elecs.mat']);
            if nanmean(maxSNRs(2,ismember(maxSNRs(1,:), elecs))) < minSNR
                continue;
            end
            
            lowBand = (rls_data(1).amp{iTrial});
            highBand = (rls_data(2).amp{iTrial});
            
            % account for trials that weren't recorded properly (mistakes)
            if isempty(lowBand)
                continue;
            end
            
            % Demean the RLS data
            lowBand = lowBand - nanmean(lowBand);
            highBand = highBand - nanmean(highBand);
            
            % Get perceptual epochs for this trial
            trialEpochs = simSchedule(simSchedule(:,1) == iTrial,:);
            
            % Plot RLS traces to visualize segments
            if strcmp(segPlotOrNot, 'yes')
                
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
            
            % partition RLS data into epochs based on sim schedule
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
                
                eMid = eStart + ((eEnd - eStart)/ 2);
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
                
                % Collect segment in raster matrix
                seg = segStart:segEnd;
                lowFreqRast = [lowFreqRast; lowBand(seg)];
                highFreqRast = [highFreqRast; highBand(seg)];
                
                % Shade segment in RLS plot and store its label
                if strcmp(segPlotOrNot, 'yes')                   
                    horiz = [seg(1)/sampRate seg(end)/sampRate];
                    vertUp = [uShadeLim uShadeLim];
                    vertDown = [lShadeLim lShadeLim];
                end
                
                if trialEpochs(iEpoch,5) > 0.5
                    if strcmp(segPlotOrNot, 'yes')
                        area(horiz,vertUp,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
                        area(horiz,vertDown,'facecolor','r','facealpha',0.15,'edgealpha',0.1);
                    end
                    labelVect = [labelVect; 'hi'];
                elseif trialEpochs(iEpoch,5) < -0.5
                    if strcmp(segPlotOrNot, 'yes')
                        area(horiz,vertUp,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
                        area(horiz,vertDown,'facecolor','b','facealpha',0.15,'edgealpha',0.1);
                    end
                    labelVect = [labelVect; 'lo'];
%                 else
%                     if strcmp(segPlotOrNot, 'yes')
%                         area(horiz,vertUp,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
%                         area(horiz,vertDown,'facecolor','y','facealpha',0.15,'edgealpha',0.1);
%                     end
                end
            end
        end
    end
end

%% Create test set (identical to training set) 
testSuffix = repmat('_t', length(labelVect), 1);
R = ones(length(labelVect),1);
testLabelVect = mat2cell([cell2mat(labelVect) testSuffix],R);
highFreqRast = [highFreqRast; highFreqRast];
lowFreqRast = [lowFreqRast; lowFreqRast];
labelVect = [labelVect; testLabelVect];

%% Save raster data
raster_data = lowFreqRast;
raster_labels.domFreq = labelVect;
raster_site_info = [];
save(['ML/raster_data/' parName '_lowFreq'], 'raster_data', 'raster_labels', 'raster_site_info');

raster_data = highFreqRast;
raster_labels.domFreq = labelVect;
raster_site_info = [];
save(['ML/raster_data/' parName '_highFreq'], 'raster_data', 'raster_labels', 'raster_site_info');
