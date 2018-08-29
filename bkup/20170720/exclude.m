% exclude.m creates lists of excluded participants and saves them as 
% variables in pre-processing/stratExcluded.mat and pre-processing/
% cumExcluded.mat.
%
% Note that this does not include individual trials that are excluded based
% on SNR
%
% Left out participants for Machine Learning Analysis have stratReason 'IQ
% matching'

% Stratus
stratExcluded = [106 112 122 124 125 129 130 131 132 133 135]; 
stratReasons = ...
    {'IQ matching'
    'IQ matching'
    'behavioral-outlier';
    'SNR';
    'EEG oscFreq-outlier';
    'behavioral-outlier';
    'SNR';
    'IQ matching';
    'IQ matching';
    'behavioral-outlier';
    'behavioral-outlier'};
save('pre-processing/stratExcluded.mat','stratExcluded','stratReasons')

% Cumulus
cumExcluded = [5 9 10 16 22 23 25 28 29]; 
cumReasons = ...
    {'SNR';
    'SNR';
    'SNR';
    'SNR';
    'SNR';
    'behavioral-outlier';
    'SNR';
    'SNR';
    'behavioral-outlier_wanderingEye'};
save('pre-processing/cumExcluded.mat','cumExcluded','cumReasons')