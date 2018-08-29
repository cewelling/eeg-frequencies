function [slowDurs, fastDurs, mixedDurs] = pressDurations(omniData,parName,analRunNum)
% -Pulls out durations of perceptual epochs identified from key press data into vectors
% slowDurs, fastDurs, and mixedDurs (separated into slow frequency, fast
% frequency, and mixed percept)
% -Plots average duration of each stimulus frequency percept (& mixed percept if applicable)  

slowReports = [];
fastReports = [];
mixedReports = [];

for i = 1:size(omniData,1)
    if omniData(i,7) == -1
        slowReports = [slowReports; omniData(i,:)];
    elseif omniData(i,7) == 1
        fastReports = [fastReports; omniData(i,:)];
    elseif omniData(i,7) == 0
        mixedReports = [mixedReports; omniData(i,:)];
    end
end


slowDurs = slowReports(:,4);
slowDur = mean(slowReports(:,4));
slowError = std(slowReports(:,4));

fastDurs = fastReports(:,4);
fastDur = mean(fastReports(:,4));
fastError = std(fastReports(:,4));

if ~isempty(mixedReports)
    mixedDurs = mixedReports(:,4);
    mixedDur = mean(mixedReports(:,4));
    mixedError = std(mixedReports(:,4));

    figure
    hold on
    bar(1:3,[slowDur fastDur mixedDur])
    Labels = {'14.166 Hz', '17 Hz', 'Mixed'};
    set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
    errorbar(1:3,[slowDur fastDur mixedDur],[slowError fastError mixedError],'.')
    title(['Average Percept Durations, ' parName ', run ' num2str(analRunNum)])
    ylabel('seconds')
    
else    
    figure
    hold on
    bar(1:2,[slowDur fastDur])
    Labels = {'14.166 Hz', '17 Hz'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    errorbar(1:2,[slowDur fastDur],[slowError fastError],'.')
    title(['Average Percept Durations, ' parName ', run ' num2str(analRunNum)])
    ylabel('seconds')
end

end

