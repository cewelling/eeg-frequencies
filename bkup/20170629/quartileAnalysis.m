function quartileAnalysis(parName, date)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

% load parameters
analysisParams

%% Load transitions
if exist(['transitions/finalTransitions/' parName '_' normType '.mat'], 'file')
    load(['transitions/finalTransitions/' parName '_' normType '.mat']);
else
    [ tLists, ~, ~, ~] = makeParPlot(parName, date);
end

f1Ratio = nan(1,2);
f2Ratio = nan(1,2);
hmLRatio = nan(1,2);

%% Get individual frequency traces, high minus low
for tType = 1:2 % low to high, high to low
    f1Traces = tLists{tType}.f1(all(~isnan(tLists{tType}.f1),2),:);
    f2Traces = tLists{tType}.f2(all(~isnan(tLists{tType}.f2),2),:);
%     f1Traces(isnan(f1Traces)) = [];
%     f2Traces(isnan(f2Traces)) = [];
    hMinusL = f2Traces - f1Traces;
    
    %% Establish amplitude limits for each transition
    lim1 = 0.40; % proportion of range
    lim2 = 0.60;
    
    f1q1 = min(f1Traces, [], 2) + range(f1Traces, 2)*lim1; f1q2 = min(f1Traces, [], 2) + range(f1Traces, 2)*lim2;
    f2q1 = min(f2Traces, [], 2) + range(f2Traces, 2)*lim1; f2q2 = min(f2Traces, [], 2) + range(f2Traces, 2)*lim2;
    hmLq1 = min(hMinusL, [], 2) + range(hMinusL, 2)*lim1; hmLq2 = min(hMinusL, [], 2) + range(hMinusL, 2)*lim2;
    
    %% Create first and 3rd quartile matrices
    pts = size(f1Traces, 2);
    f1q1 = repmat(f1q1, 1, pts); f1q2 = repmat(f1q2, 1, pts);
    f2q1 = repmat(f2q1, 1, pts); f2q2 = repmat(f2q2, 1, pts);
    hmLq1 = repmat(hmLq1, 1, pts); hmLq2 = repmat(hmLq2, 1, pts);
    
    %% Store proportion of time spent in interquartile range
    f1Ratio(tType) = nansum(nansum(f1Traces > f1q1 & f1Traces < f1q2))/numel(f1Traces);
    f2Ratio(tType) = nansum(nansum(f2Traces > f2q1 & f2Traces < f2q2))/numel(f2Traces);
    hmLRatio(tType) = nansum(nansum(hMinusL > hmLq1 & hMinusL < hmLq2))/numel(hMinusL);
end

%% Save average proportion for participant
meanf1Ratio = nanmean(f1Ratio);
meanf2Ratio = nanmean(f2Ratio);
meanHmLRatio = nanmean(hmLRatio);
save([indicesDir 'interQprop/' parName], 'meanf1Ratio', 'meanf2Ratio', 'meanHmLRatio');

end