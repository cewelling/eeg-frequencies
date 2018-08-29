% % Find and store transitions between percepts that meet criteria
%     
% lowToHigh = []; l2hTrials = [];
% highToLow = []; h2lTrials = [];
% mixedToHigh = []; m2hTrials = [];
% mixedToLow = []; m2lTrials = [];
% highToMixed = []; h2mTrials = [];
% lowToMixed = []; l2mTrials = [];
% 
% for iEpoch = 1:(size(epochs,1) - 1) %don't check the transition out of the last epoch
%     thisEpoch = epochs(iEpoch, :);
%     nextEpoch = epochs(iEpoch + 1, :);
%     if nextEpoch(2) - thisEpoch(3) <= gapMax && nextEpoch(1) == thisEpoch(1) % epochs must be in the same trial
%         if thisEpoch(5) == -1
%             if thisEpoch(3) - thisEpoch(2) >= domMin
%                 if nextEpoch(5) == 1
%                     if nextEpoch(3) - nextEpoch(2) >= domMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         lowtoHigh = [lowToHigh; transPt];
%                         l2hTrials = [l2hTrials; thisEpoch(1)];
%                     end
%                 elseif nextEpoch(5) == 0
%                     if nextEpoch(3) - nextEpoch(2) >= mixedMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         lowToMixed = [lowToMixed; transPt];
%                         l2mTrials = [l2mTrials; thisEpoch(1)];
%                     end
%                 end
%             end
%             
%         elseif thisEpoch(5) == 1
%             if thisEpoch(3) - thisEpoch(2) >= domMin
%                 if nextEpoch(5) == -1
%                     if nextEpoch(3) - nextEpoch(2) >= domMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         highToLow = [highToLow; transPt];
%                         h2lTrials = [h2lTrials; thisEpoch(1)];
%                     end
%                 elseif nextEpoch(5) == 0
%                     if nextEpoch(3) - nextEpoch(2) >= mixedMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         highToMixed = [highToMixed; transPt];
%                         h2mTrials = [h2mTrials; thisEpoch(1)];
%                     end
%                 end
%             end
%             
%         elseif thisEpoch(5) == 0
%             if thisEpoch(3) - thisEpoch(2) >= mixedMin
%                 if nextEpoch(5) == -1
%                     if nextEpoch(3) - nextEpoch(2) >= domMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         mixedToLow = [mixedToLow; transPt];
%                         m2lTrials = [m2lTrials; thisEpoch(1)];
%                     end
%                 elseif nextEpoch(5) == 1
%                     if nextEpoch(3) - nextEpoch(2) >= domMin
%                         transPt = mean([thisEpoch(3) nextEpoch(2)]);
%                         mixedToHigh = [mixedToHigh; transPt];
%                         m2hTrials = [m2hTrials; thisEpoch(1)];
%                     end
%                 end
%             end
%         end
%     end
% end


% AlinaSNRsFreq1 = [];
% AlinaSNRsFreq2 = [];
% for parNum = 133:136
%     if exist('indices/SNR/stratus' parNum '_run1_oz.mat')
%     load(['indices/SNR/stratus' parNum '_run1_oz.mat')
%     newSNRfreq1 = meanSNRs(1).value;
%     newSNRfreq2 = meanSNRs(2).value;
%     AlinaSNRsFreq1 = [AlinaSNRsFreq1 newSNRfreq1];
%     AlinaSNRsFreq2 = [AlinaSNRsFreq2 newSNRfreq2];
% end
% 
% AlinaSNRmeanFreq1 = mean(AlinaSNRsFreq1);
% AlinaSNRmeanFreq2 = mean(
    
% load & gather SNRs
