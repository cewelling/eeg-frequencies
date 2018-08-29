clearvars

%load parameters
analysisParams

figure

numFormat = '%02d'; % for participant number

% space to store individual's characteristic freq of oscillation
stratusFreqs = [];
cumFreqs = [];

for groupNum = [0 1]
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
    
    for iPar = 1:length(parNums)
        parName = [group num2str(parNums(iPar), numFormat)];
        try
            [freq, meanfft] = rateFFT(parName);
        catch
            continue;
        end
        ffts = [ffts; meanfft];
        [~,I] = max(meanfft);
        if groupNum == 0 
            stratusFreqs = [stratusFreqs freq(I)];
        else
            cumFreqs = [cumFreqs freq(I)];
        end
    end
    
    meantrace = mean(ffts);
    errortrace = ste(ffts);
    
    % plot average fft trace
    mseb(freq,meantrace,errortrace,lineProps,1);
    %plot(freq, meantrace)
    [~,I] = max(meantrace);
    vline(freq(I));
    hold on
    
end

% fft trace plot formatting
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
