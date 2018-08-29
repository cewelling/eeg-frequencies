% exclude.m creates lists of excluded participants and saves them as 
% variables in pre-processing/stratExcluded.mat and pre-processing/
% cumExcluded.mat.
%
% Note that this does not include individual trials that are excluded based
% on SNR
%
% This only includes behavioral outliers. Oscillation frequency outliers
% are excluded in indexAnal2.m. 

% Stratus
stratExcluded = [106 112 122 129 131 132 133 134 135]; %[122 125 129 133 135]; %[125 133]; %[116 122 125 129 133 135];
stratReasons = ...
    {'IQ matching'
    'IQ matching'
    'behavioral-outlier';
    'behavioral-outlier';
    'IQ matching';
    'IQ matching';
    'behavioral-outlier';
    'no IQ score';
    'behavioral-outlier'};
save('pre-processing/stratExcluded.mat','stratExcluded','stratReasons')

% Cumulus
cumExcluded = [23 29]; %[6 16 29]; %[3 6 14 16 18 21 25 29]; 
cumReasons = ...
    {'behavioral-outlier';
    'behavioral-outlier_wanderingEye'};
save('pre-processing/cumExcluded.mat','cumExcluded','cumReasons')