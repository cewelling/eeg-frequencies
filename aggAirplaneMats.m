% After running airplane mats for a bunch of latencies, aggregate them into
% a single 3D matrix

airplaneDir = 'ML/airplaneMats/';

% Establish which runs to aggregate
if ~exist('firstRun', 'var')
    firstRun = 4;
    lastRun = 8;
end

fullMat = [];
allRivLats = [];
for thisRun = firstRun:lastRun
    runID = sprintf('%03u',thisRun);
    files = dir([airplaneDir 'run' runID '*.mat']);
    if ~isempty(files)
        load([airplaneDir files(1).name]);
        fullMat = cat(1, fullMat, airplaneMat);
        allRivLats = [allRivLats rivLats];
%     else
%         fullMat = cat(1, fullMat, nan(1, size(fullMat, 2), size(fullMat,3)));
%         allRivLat = [allRivLats NaN];
     end
end

% for consistency with matrices that weren't aggregated
airplaneMat = fullMat;
rivLats = allRivLats; 

% Now we have a 3D matrix with:
% Dim 1: latency
% Dim 2: training seg
% Dim 3: testing timepoint

save(['ML/airplaneMats/runs' num2str(firstRun) '-' num2str(lastRun)], 'airplaneMat', 'rivLats', 'segTime', 'teTime');

