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
    parNums = groupParNums{iGroup};
    group = groupCodes{iGroup};
    
    % line colors for plotting
    if iGroup == 1
        lineProps.col = {'b'};
    else
        lineProps.col = {'r'};
    end
    
    ffts = [];
    freq = [];
    for parNum = parNums
        parName = [group num2str(parNum, numFormat)];
        
        % Get EEG file info to determine whether participant exists
        EEGfiles = dir([eegDir parName '*']);
        % ...if the participant exists...
        if isempty(EEGfiles)
            continue;
        end
        
        % get oscillation frequency, mean trace for participant
        [freq, meanfft] = rateFFT(parName);
        
        % get frequency axis (should be the same for each participant)
        if ~exist('freqAxis', 'var') && ~isempty(freq)
            freqAxis = freq;
        end
        
        % collect fft traces
        ffts = [ffts; meanfft];
        
        % collect oscillation frequencies (stored in rateFFT.m)
        load([indicesDir 'rateFreqs/' parName '.mat']);
        if iGroup == 1 
            stratusFreqs = [stratusFreqs oscFreq];
        else
            cumFreqs = [cumFreqs oscFreq];
        end
    end
    
    meantrace{iGroup} = mean(ffts);
    errortrace{iGroup} = ste(ffts);    
end

% plot average fft trace
figure
for iGroup = 1:length(groupParNums)
    % plot average trace
    mseb(freqAxis,meantrace{iGroup},errortrace{iGroup},lineProps,1);
    
    % calculate CDF
    CDF = cumsum(meantrace{iGroup});
    CDF = CDF/max(CDF);
    
    % find the index of the half max
    [~, hMaxI] = min(abs(.5 - CDF));
    
    % Mark the frequency of the half max
    if iGroup == 1
        vline(freqAxis(hMaxI), 'b');
    else
        vline(freqAxis(hMaxI), 'r');
    end
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
