clearvars 

channel = 29;
runTypes = {'dartRival', 'dartSim'};

for iRunType = [1 2];
    runType = runTypes{iRunType};
    traceFiles = dir(['*' num2str(channel) '_' runType '.mat']);
    lowTraces = [];
    highTraces = [];
    labelTraces = [];
    for i = 1:length(traceFiles)
        load(traceFiles(i).name);
        lowTraces = [lowTraces; nanmean(lowMat)];
        highTraces = [highTraces; nanmean(highMat)];
        labelTraces = [labelTraces; nanmean(labelMat)];
    end
    
    % plot traces
    figure;
    plot(nanmean(lowTraces));
    hold on;
    plot(nanmean(highTraces));
    plot(nanmean(labelTraces));
    title(['stratus ' runType ': high to low'])
    vline(5122/2);
end