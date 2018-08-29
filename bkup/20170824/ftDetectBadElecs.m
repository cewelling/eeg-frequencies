function [ badElecs ] = ftDetectBadElecs(checkNoise_data, runName)
% ftDetectBadElecs facilitates visual detection of noisy electrodes prior 
% to fft analysis.

%% Create histograms for noisy electrode detection
kurMatrix = [];
meanMatrix = [];
varMatrix = [];
for iTrial = 1:length(checkNoise_data.trial)
    trialData = checkNoise_data.trial{iTrial};
    %for elec = 1:size(trialData,1)
    for elec = 1:32
        kur = kurtosis(trialData(elec,:));
        mn = mean(trialData(elec,:));
        vr = var(trialData(elec,:));
        kurMatrix = [kurMatrix kur];
        meanMatrix = [meanMatrix mn];
        varMatrix = [varMatrix vr];
    end
end

figure
hist(kurMatrix, 250);
title([runName ' kurtosis'])
drawnow
adjustHist(kurMatrix)
fprintf('choose min \n');
[kurMin,~] = ginput(1);
fprintf('kurtosis min selected \n');
fprintf('choose max \n');
[kurMax,~] = ginput(1);
fprintf('kurtosis max selected \n');
close gcf

figure
hist(meanMatrix, 250)
title([runName ' mean'])
drawnow
adjustHist(meanMatrix)
fprintf('choose min \n');
[meanMin,~] = ginput(1);
fprintf('mean min selected \n');
fprintf('choose max \n');
[meanMax,~] = ginput(1);
fprintf('mean max selected \n');
close gcf

figure
hist(varMatrix, 250)
title([runName ' variance'])
drawnow
adjustHist(varMatrix)
fprintf('choose min \n');
[varMin,~] = ginput(1);
fprintf('variance min selected \n');
fprintf('choose max \n');
[varMax,~] = ginput(1);
fprintf('variance max selected \n');
close gcf

%% Find noisy electrodes
badElecs = {};

for iTrial = 1:length(checkNoise_data.trial)
    trialData = checkNoise_data.trial{iTrial};
    badCount = 1;
    %         for iElec = 1:size(trialData,1)
    for iElec = 1:32
        kur = kurtosis(trialData(iElec,:));
        vr = var(trialData(iElec,:));
        mn = mean(trialData(iElec,:));
        if kur < kurMin || kur > kurMax ...
                || vr < varMin || vr > varMax ...
                || mn < meanMin || mn > meanMax
            badElecs{iTrial}(badCount).num = iElec;
            badElecs{iTrial}(badCount).kur = kur;
            badElecs{iTrial}(badCount).var = vr;
            badElecs{iTrial}(badCount).mean = mn;
            badCount = badCount + 1;
        end
    end
end

end

