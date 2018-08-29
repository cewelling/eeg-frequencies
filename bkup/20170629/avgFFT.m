clearvars -except fftCons optFFT findCorr ampClust

%load parameters
analysisParams
groupAnalysisParams

numFormat = '%02d'; % for participant number

% space to store individual's characteristic freq of oscillation
stratusFreqs = [];
cumFreqs = [];

% space to store group traces
meanTraces = {};
errorTraces = {};

for iGroup = 1:length(groupParNums)
    if groupNum == 0
        parNums = [105:140]; %[105:140]; % ex/ Stratus100's parNum is 100
        group = 'stratus'; % 'stratus' or 'cumulus'
        lineProps.col = {'b'};
    else
        parNums = [1:27]; %[1:27]; % ex/ Stratus100's parNum is 100
        group = 'cumulus'; % 'stratus' or 'cumulus'
        lineProps.col = {'r'};
    end
    ffts = [];
    freq = [];
    for parNum = groupParNums{iGroup}
        parName = [group num2str(parNums(iPar), numFormat)];
        
        % Limit to participants with at least 3 switches (for fft analysis)
        if exist([bResultsDir parName '_darts_dartRival_Switches.txt'])
            switchNum = load([bResultsDir parName '_darts_dartRival_Switches.txt']);
            % Limit number of switches (for purpose of fft analysis)
            if switchNum < 3
                continue;
            end
        end
        
        % get oscillation frequency, mean trace for participant
        [freq, meanfft] = rateFFT(parName);
        ffts = [ffts; meanfft];
        [~,I] = max(meanfft);
        if groupNum == 0 
            stratusFreqs = [stratusFreqs freq(I)];
        else
            cumFreqs = [cumFreqs freq(I)];
        end
    end
    
    meantrace{groupNum + 1} = mean(ffts);
    errortrace{groupNum + 1} = ste(ffts);    
end

% plot average fft trace
figure
for groupNum = [0 1]
mseb(freq,meantrace{groupNum + 1},errortrace{groupNum+1},lineProps,1);
plot(freq, meantrace{groupNum + 1})
[~,I] = max(meantrace{groupNum + 1});
vline(freq(I));
hold on
end
xlabel('frequency (hz)')
ylabel('power')
xlim([0 2])
legend('stratus', 'cumulus')
title('Average FFT Traces')

% bar graph
figure
h = bar([1 2], [nanmean(stratusFreqs) nanmean(cumFreqs)]);
hold on
scatterX = [ones(1, length(stratusFreqs)) (ones(1, length(cumFreqs)) + 1)];
scatter(scatterX, [stratusFreqs cumFreqs],'filled')
errorbar([1 2], [nanmean(stratusFreqs) nanmean(cumFreqs)], [ste(stratusFreqs) ste(cumFreqs)], 'k.','LineWidth',2);
xlabel('Controls                                                                   ASC  ')
ylabel('Oscillation frequency (Hz)')

% difference between controls and ASCs?
[H,P] = ttest2(stratusFreqs,cumFreqs);
text(1.5, 0.2,['p = ' num2str(P)])
