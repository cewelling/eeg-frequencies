% Plot airplane matrix
% x-axis: training segment (time before end of epoch (s))
% y-axis: rivalry latency (s)
% Color: 
%    - Matrix 1: accuracy of peak to the left of 0
%    - Matrix 2: accuracy of peak to the right of 0
%    - Matrix 3: average of both peaks
%
% Dependencies: aggAirplaneMats.m
% Called from: 

clearvars 

% airplane matrix run(s) to plot
firstRun = 38;
lastRun = 38;

% load airplane matrix
if firstRun ~= lastRun
    if ~exist(['ML/airplaneMats/runs' num2str(firstRun) '-' num2str(lastRun) '.mat'], 'file')
        aggAirplaneMats;
    end
    load(['ML/airplaneMats/runs' num2str(firstRun) '-' num2str(lastRun) '.mat'])
else
    runID = sprintf('%03u', firstRun);
    airplaneMatFile = dir(['ML/airplaneMats/run' runID '*.mat']);
    if length(airplaneMatFile) ~= 1
        error(['More than one file named run ' runID]);
    end
    load(['ML/airplaneMats/' airplaneMatFile(1).name])
end

% Remember, the dimensions of the airplane matrix are:
% Dim 1: latency                
% Dim 2: training seg
% Dim 3: testing timepoint

% plotting
figure;

clims = [0.5 1]; % color limits

% accuracy of left peak
subplot(1,3,1)
leftHalf = 1:floor(size(airplaneMat, 3)/2);
leftMat = squeeze(nanmax(airplaneMat(:,:,leftHalf), [], 3));
imagesc(segTime, rivLats, leftMat, clims);
title('Left Peak Accuracies')
%xlabel('Training segment - (s) before epoch end')
ylabel('Rivalry latency (s)')
colormap(jet);
colorbar
set(gca,'fontsize',14)

% accuracy of right peak
subplot(1,3,2)
rightHalf = (ceil(size(airplaneMat, 3)/2) + 1):size(airplaneMat, 3);
rightMat = squeeze(nanmax(airplaneMat(:,:, rightHalf), [], 3));
imagesc(segTime, rivLats, rightMat, clims);
title('Right Peak Accuracies')
xlabel('Training segment - (s) before epoch end')
ylabel('Rivalry latency (s)')
colormap(jet);
colorbar
set(gca,'fontsize',14)

% accuracy of average
subplot(1,3,3)
leftAndRight = cat(3, leftMat, rightMat);
avgMat = squeeze(nanmean(leftAndRight, 3));
imagesc(segTime, rivLats, avgMat, clims);
title('Average of Right and Left Peaks');
%xlabel('Training segment - (s) before epoch end')
ylabel('Rivalry latency (s)')
colormap(jet);
colorbar
set(gca,'fontsize',14)

