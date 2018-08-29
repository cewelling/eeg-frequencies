clear all; close all; clc;
% setFigProps

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runsOfInt = [1 2 3]; %[7 8 11 12]; %[1 4 5 6];
par = [105];

lowNormFreqVals = [];
highNormFreqVals = [];

winLengthRLS = 1000; %time in ms

loopCount = 0;
% for par = pars
    
%     if par == 100
%         runsOfInt = [2 10];
%     else
%         runsOfInt = [4 5 6];
%     end

for analRunNum = runsOfInt %[7 8 11 12] % [1 4 5 6]

%%%%% Set up: Things to change %%%%%%%%%%%%%%%
% analRunNum = 6;
parName = ['stratus' num2str(par)];
% date = '4_15';
% fileNames = dir(['../../runScripts/Results/' parName '*' date '*']);
fileNames = dir(['../../runScripts/Results/' parName '*']);
legend = importdata(['pilotlegends/' parName '.txt']); %pilot-' date '-caroline.txt
stimulation_frequencies = [14.16 17]; %[17 21.25] Note: this is redefined later....
intermodulation_frequencies = 11.3;

%%%%% Set up: Always done %%%%%%%%%%%%%%%%%%%%%
fileNames = {fileNames.name}; filenames{1}(1).name = char(fileNames(analRunNum));
stimulation_frequencies = legend.data(analRunNum,1:2);
numTrials = legend.data(analRunNum,3);
trialDur = legend.data(analRunNum,4); 
defaultAnalParams %Load default anal params

% stimulation_frequencies = [7.0833 8.5000];

stimulation_frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omniData = runKeyAnal(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,minEpoch);

multiFreqSeries = [];
trialTime = [];

[multiFreqSeries, trialTime] = compareRLSfast(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, numTrials, trialDur,parName,legend, omniData,winLengthRLS);


%Normalize the frequency band time series
normalFreqSeries = [];

for iTrial = 1:size(multiFreqSeries,2)
    
    for iFreq = 1:size(multiFreqSeries{1,iTrial},2);
        
        RLSamp = cell2mat(multiFreqSeries{1,iTrial}(1,iFreq));
        
        normRLSamp = 2*(RLSamp - min(RLSamp))/(max(RLSamp) - min(RLSamp)) - 1;
        
%         normRLSamp2 = (RLSamp - ((max(RLSamp)+min(RLSamp))/2))/((max(RLSamp)-min(RLSamp))/2);
%         normRLSamp = zscore(RLSamp);
        
        normalFreqSeries{1,iTrial}{1,iFreq} = normRLSamp;
        
    end
    
end

for iTrial = 1:size(normalFreqSeries,2)
    for iFreq = 1:size(normalFreqSeries{1,iTrial},2)
        normalFreqVals = cell2mat(normalFreqSeries{1,iTrial}(1,iFreq));
        if iFreq == 1
            lowNormFreqVals = [lowNormFreqVals; normalFreqVals];
        elseif iFreq == 2
            highNormFreqVals = [highNormFreqVals; normalFreqVals];
        end
    end
end
        


loopCount = loopCount+1;
close all
end
% end

figure
scatter(lowNormFreqVals,highNormFreqVals,'k','.')
axis([-1.5 1.5 -1.5 1.5])

%downsample
sparseLowNorm = downsample(lowNormFreqVals,50);
sparseHighNorm = downsample(highNormFreqVals,50);
r = corr2(sparseLowNorm,sparseHighNorm);
r = round(r,3);

figure
scatter(sparseLowNorm,sparseHighNorm,'k','.')
hold on
title(parName)
axis([-1.5 1.5 -1.5 1.5])
h = lsline;
text(0.75, 1, ['r = ' num2str(r)],'FontSize',14)

set(h,'color','r')
set(h,'LineWidth',1)
