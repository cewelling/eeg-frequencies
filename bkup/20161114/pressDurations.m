function [slowDurs, fastDurs, mixedDurs] = pressDurations(omniData,parName,analRunNum)
% -Pulls out durations of perceptual epochs identified from key press data into vectors
% slowDurs, fastDurs, and mixedDurs (separated into slow frequency, fast
% frequency, and mixed percept)
% -Plots average duration of each stimulus frequency percept (& mixed percept if applicable)  

slowReports = [];
fastReports = [];
mixedReports = [];

for i = 1:size(omniData,1)
    if omniData(i,5) == -1
        slowReports = [slowReports; omniData(i,:)];
    elseif omniData(i,7) == 1
        fastReports = [fastReports; omniData(i,:)];
    elseif omniData(i,7) == 0
        mixedReports = [mixedReports; omniData(i,:)];
    end
end

% Calculate and plot mean duration and error for each percept
numBars = 0;
barContent = [];
barLabels = {};
errorBars = [];

if ~isempty(slowReports)
    slowDurs = slowReports(:,4);
    slowDur = mean(slowReports(:,4));
    slowError = std(slowReports(:,4));
    numBars = numBars + 1;
    barContent = [barContent slowDur];
    barLabels = [barLabels, '5.67 Hz'];
    errorBars = [errorBars slowError];
else
    slowDurs = NaN;
end
    
if ~isempty(fastReports)
    fastDurs = fastReports(:,4);
    fastDur = mean(fastReports(:,4));
    fastError = std(fastReports(:,4));
    numBars = numBars + 1;
    barContent = [barContent fastDur];
    barLabels = [barLabels, '8.5 Hz'];
    errorBars = [errorBars fastError];
else
    fastDurs = NaN;
end

if ~isempty(mixedReports)
    mixedDurs = mixedReports(:,4);
    mixedDur = mean(mixedReports(:,4));
    mixedError = std(mixedReports(:,4));
    numBars = numBars + 1;
    barContent = [barContent mixedDur];
    barLabels = [barLabels, 'Mixed'];
    errorBars = [errorBars mixedError];
else 
    mixedDurs = NaN;
end

    figure
    hold on
    bar(1:numBars, barContent)
    Labels = barLabels;
    set(gca, 'XTick', 1:numBars, 'XTickLabel', Labels);
    errorbar(1:numBars, barContent, errorBars,'.')
    title(['Average Percept Durations, ' parName ', run ' num2str(analRunNum)])
    ylabel('seconds')

end

