% aggregate accMats (which have accuracies for all participants, but just
% one segment and 1 latency) into airplaneMats (which have accuracies
% averaged across participants for all segments, but only one latency.

%% Edit this section with desired settings
group = 'stratus';

trainSet = 'sim';
testSet = 'rivalTransitions';

% latencyNum to aggregate
latNum = 1;
rivLats = 0; % s
rivLat = rivLats; % since we are saving one file per latency

% segments tested
segTime = 0.25:0.05:0.75;

% runID for filename
runID = '002';

%% Aggregate files into airplane matrix
files = dir([group '_train-' trainSet '_test-' testSet '_lat' num2str(latNum) '_seg*']); 

for i = 1:length(files)
    load([group '_train-' trainSet '_test-' testSet '_lat' num2str(latNum) '_seg' num2str(i) '.mat'])
    if ~exist('airplaneMat', 'var')
        airplaneMat = nan(1,length(segTime),length(teTime));
    end
    airplaneMat(1,i,:) = nanmean(accuracyMat);
end

%% Save file

save(['../airplaneMats/run' runID '_' group '_' num2str(size(accuracyMat, 1)) 'pars_train-' trainSet '_test-' testSet '_' num2str(round(1000*rivLat)) 'msLat'], 'airplaneMat', 'rivLats', 'segTime','teTime')

    