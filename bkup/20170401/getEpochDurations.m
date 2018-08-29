function getEpochDurations(parName, runName, date)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

% calculate epoch durations


% load parameters
analysisParams

% get list of participant's perceptual epochs
[epochs, simSchedule] = getEpochs(parName, runName, date);

slowReports = [];
fastReports = [];
mixedReports = [];

for i = 1:size(epochs,1)
    if epochs(i,5) == -1
        slowReports = [slowReports; epochs(i,:)];
    elseif epochs(i,5) == 1
        fastReports = [fastReports; epochs(i,:)];
    elseif epochs(i,5) == 0
        mixedReports = [mixedReports; epochs(i,:)];
    end
end

% Calculate and plot mean duration and error across percepts (one plot has
% one run)
numBars = 0;
barContent = [];
barLabels = {};
errorBars = [];

if ~isempty(slowReports)
    slowDurs = slowReports(:,3) - slowReports(:,2);
    slowDur = mean(slowDurs);
    slowError = std(slowDurs);
    numBars = numBars + 1;
    barContent = [barContent slowDur];
    barLabels = [barLabels, '5.67 Hz'];
    errorBars = [errorBars slowError];
else
    slowDurs = NaN;
    slowDur = NaN;
    slowError = NaN;
end
    
if ~isempty(fastReports)
    fastDurs = fastReports(:,3) - fastReports(:,2);
    fastDur = mean(fastDurs);
    fastError = std(fastDurs);
    numBars = numBars + 1;
    barContent = [barContent fastDur];
    barLabels = [barLabels, '8.5 Hz'];
    errorBars = [errorBars fastError];
else
    fastDurs = NaN;
    fastDur = NaN;
    fastError = NaN;
end

if ~isempty(mixedReports)
    mixedDurs = mixedReports(:,3) - mixedReports(:,2);
    mixedDur = mean(mixedDurs);
    mixedError = std(mixedDurs);
    numBars = numBars + 1;
    barContent = [barContent mixedDur];
    barLabels = [barLabels, 'Mixed'];
    errorBars = [errorBars mixedError];
else
    mixedDurs = NaN;
    mixedDur = NaN;
    mixedError = NaN;
end

figure
hold on
bar(1:numBars, barContent)
Labels = barLabels;
set(gca, 'XTick', 1:numBars, 'XTickLabel', Labels);
errorbar(1:numBars, barContent, errorBars,'.')
title(['Average Percept Durations, ' parName ',' runName])
ylabel('seconds')

% Save durations, supression indices
durData = [slowDur slowError; fastDur fastError; mixedDur mixedError];
save ([indicesDir 'EpochDurations/' parName '_' runName], 'durData');
domSum = nansum(slowDurs) + nansum(fastDurs);
suppIndex = domSum / (nansum(mixedDurs) + domSum);
save ([indicesDir 'suppIndices/' parName '_' runName], 'suppIndex');

end