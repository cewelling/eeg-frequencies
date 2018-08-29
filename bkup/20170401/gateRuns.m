function [ runsIncluded ] = gateRuns( input_args )
%UNTITLED19 Summary of this function goes here
%   Detailed explanation goes here

runList = [];
for currentRun = runIndices
    if size(Legend.data,1) < currentRun
        continue;
    end
    % Get the appropriate EEG file for this run
    EEGfile = [results_directory EEGfiles(currentRun).name];
    SNRs = runFFT(EEGfile, stimFreqs, currentRun, numTrials, trialDur, parName, iPar);
    if length(SNRs) > 2
        error('Choose just one noise window, electrode set for %s run %d SNR', parName, currentRun)
    end
    if mean([SNRs.value]) > 3
        runList = [runList currentRun];
    end
end

end

