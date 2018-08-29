function [transitions] = compareEpochs(omniData, multiFreqSeries, trialTime, analRunNum,epochHalf,normType)
 
%normType = 'z' or 'norm' or 'none'

% Transition codes to look for in the omniData, and their labels
caseCode = [-1 1; 1 -1; 0 1; 0 -1; 1 0; -1 0];
caseType = {'low to high'; 'high to low'; 'mixed to high'; 'mixed to low'; 'high to mixed'; 'low to mixed'};

% Set the time requirements for consideration in the transition analysis
dominantMin = 1.25;
mixedMin = 1;
betweenReportMax = 0.25; % maximum time allowed between button-reported states

% epochHalf = 1.5; % this is chosen in the transitionAnal code now

% Preallocate space
caseTimes = zeros(50,6);
caseTrials = zeros(50,6);


% loop through the button report data to find transitions between percepts
% which meet criteria
    % more specifically, for each type of possible transition, set the
    % criteria and search through the omniData for examples of these
    % transitions

for iCase = 1:size(caseCode,1)
    
    % for the transition case, set the minimum durations
    
    if caseCode(iCase,1) ~= 0
        firstStateMin = dominantMin;
    else
        firstStateMin = mixedMin;
    end
    
    if caseCode(iCase,2) ~= 0
        secondStateMin = dominantMin;
    else
        secondStateMin = mixedMin;
    end
    
    caseEpochs = []; %preallocate for accruing center points
    caseEpochTrials = [];
    
    for iEpoch = 1:(size(omniData,1) - 1) %don't check the transition out of the last epoch

        % for each line in omniData, check to see if epoch meets criteria:
        % (column 7 contains "reporting fast, slow, or mixed percept")
        if omniData(iEpoch, 7) == caseCode(iCase,1) && ...  %the first of the two states is correct
                omniData(iEpoch, 4) > firstStateMin && ...  %the first of the two states is sufficiently long
                omniData((iEpoch+1),7) == caseCode(iCase,2) && ...  %the second state is correct
                omniData((iEpoch+1), 4) > secondStateMin && ... %the second state is sufficiently long
                0 < (omniData((iEpoch+1),2) - omniData(iEpoch,3)) && ...
                (omniData((iEpoch+1),2) - omniData(iEpoch,3)) < betweenReportMax %the time between reports is sufficiently small, and not negative
            
            centerPt = mean([omniData(iEpoch,3) omniData((iEpoch+1),2)]); % if the epoch transition meets criteria, store the center time point
            caseEpochs = [caseEpochs; centerPt];
            caseEpochTrials = [caseEpochTrials; omniData(iEpoch,1)];
            
        end
        
    end
    
    caseTimes(1:length(caseEpochs),iCase) = caseEpochs;
    caseTrials(1:length(caseEpochs),iCase) = caseEpochTrials;
    
end


%Normalize the RLS data

normalFreqSeries = [];

for iTrial = 1:size(multiFreqSeries,2)
    
    for iFreq = 1:size(multiFreqSeries{1,iTrial},2);
        
        RLSamp = cell2mat(multiFreqSeries{1,iTrial}(1,iFreq));
        
        % take normalizing/standardizing step
        if strcmp(normType, 'norm')
            normRLSamp = 2*(RLSamp - min(RLSamp))/(max(RLSamp) - min(RLSamp)) - 1;
        elseif strcmp(normType, 'z')
            normRLSamp = zscore(RLSamp);
        elseif strcmp(normType, 'none')
            normRLSamp = RLSamp;
        end
        
        normalFreqSeries{1,iTrial}{1,iFreq} = normRLSamp;
        
    end
    
end

% Get the RLS time-amplitude data for the times corresponding to the
% transitions

halfLength = 512*epochHalf; % Remember we have 512 data points per second
fullLength = halfLength*2;

for iCase = 1:size(caseTrials,2) %for each type of transition
    
    thisCaseTrials = nonzeros(caseTrials(:,iCase)); %list of trials containing this transition
    
    f1Epochs = [];
    f2Epochs = [];
    
    for iEpoch = 1:length(thisCaseTrials)
        
        transitionTime = caseTimes(iEpoch,iCase);
        
        if ((transitionTime - epochHalf) >= 0) && ((transitionTime + epochHalf) <= 30) % epoch isn't cut off by beginning or end of trial
        
            RLSf1 = cell2mat(normalFreqSeries{1,thisCaseTrials(iEpoch)}(1,1));
            epochRLSf1 = RLSf1(trialTime>=transitionTime-epochHalf & trialTime<=transitionTime+epochHalf)';
            
            RLSf2 = cell2mat(normalFreqSeries{1,thisCaseTrials(iEpoch)}(1,2));
            epochRLSf2 = RLSf2(trialTime>=transitionTime-epochHalf & trialTime<=transitionTime+epochHalf)';
            
        elseif (transitionTime - epochHalf) < 0 % epoch is cut off by beginning, pad beginning with NaNs
            
            RLSf1 = cell2mat(normalFreqSeries{1,thisCaseTrials(iEpoch)}(1,1));
            epochRLSf1 = RLSf1(trialTime>=transitionTime-epochHalf & trialTime<=transitionTime+epochHalf)';

            RLSf2 = cell2mat(normalFreqSeries{1,thisCaseTrials(iEpoch)}(1,2));
            epochRLSf2 = RLSf2(trialTime>=transitionTime-epochHalf & trialTime<=transitionTime+epochHalf)';

            padSize = fullLength - length(epochRLSf1);

            epochRLSf1 = [zeros(1,padSize), epochRLSf1];
            epochRLSf1(epochRLSf1 == 0) = NaN;

            epochRLSf2 = [zeros(1,padSize), epochRLSf2];
            epochRLSf2(epochRLSf2 == 0) = NaN;
            
        elseif (transitionTime + epochHalf) > 30 % epoch is cut off by end, pad end with NaNs
                
            RLSf1 = cell2mat(normalFreqSeries{1,thisCaseTrials(iEpoch)}(1,1));
            epochRLSf1 = RLSf1(trialTime>=transitionTime-epochHalf & trialTime<=transitionTime+epochHalf)';

            RLSf2 = cell2mat(normalFreqSeries{1,thisCaseTrials(iEpoch)}(1,2));
            epochRLSf2 = RLSf2(trialTime>=transitionTime-epochHalf & trialTime<=transitionTime+epochHalf)';

            padSize = fullLength - length(epochRLSf1);

            epochRLSf1 = [epochRLSf1,zeros(1,padSize)];
            epochRLSf1(epochRLSf1 == 0) = NaN;

            epochRLSf2 = [epochRLSf2,zeros(1,padSize)];
            epochRLSf2(epochRLSf2 == 0) = NaN;
        end
        
        
            try
            
            f1Epochs = [f1Epochs; epochRLSf1];
            f2Epochs = [f2Epochs; epochRLSf2];

        catch % handle concat errors (duration calculations do not always line up with data points)
            if length(epochRLSf1) > size(f1Epochs,2)
                epochRLSf1 = epochRLSf1(1:size(f1Epochs,2));
                epochRLSf2 = epochRLSf2(1:size(f1Epochs,2));
            elseif length(epochRLSf1) < size(f1Epochs,2)
                epochRLSf1 = [epochRLSf1, epochRLSf1(end)*ones(1,(size(f1Epochs,2)-length(epochRLSf1)))];
                epochRLSf2 = [epochRLSf2, epochRLSf2(end)*ones(1,(size(f1Epochs,2)-length(epochRLSf2)))];
            end
            
            f1Epochs = [f1Epochs; epochRLSf1];
            f2Epochs = [f2Epochs; epochRLSf2];
            
        end
        
    end
    
    transitions(iCase).type = caseType(iCase);
    transitions(iCase).f1 = f1Epochs;
    transitions(iCase).f2 = f2Epochs;
        
end

% Plotting
for iCase = 1:size(transitions,2) % for each type of transition

    if isempty(transitions(iCase).f1)
        disp(['no instances of' caseType(iCase)])

    elseif size(nonzeros(caseTrials(:,iCase)),1) < 2
        disp(['not enough instances of' caseType(iCase)])

    else
        meanF1trace = mean(transitions(iCase).f1,1);
        errorF1trace = ste(transitions(iCase).f1);

        meanF2trace = mean(transitions(iCase).f2,1);
        errorF2trace = ste(transitions(iCase).f2);

        figure
        title(caseType(iCase))
        hold on
        mseb([(-epochHalf+(1/512)):(1/512):epochHalf],[meanF1trace; meanF2trace],[errorF1trace; errorF2trace],[],1);
        vline(0,'k','transition')
        xlabel('Time from button press (s)');
        ylabel('Normalized amplitude');

    end

end

end