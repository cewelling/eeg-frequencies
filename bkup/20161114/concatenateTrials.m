% Concatonate artifact-free trials
function concatData = concatenateTrials(data)

dataConcat = [];
timeConcat = [];

for i = 1:size(data.trial,2)
    dataConcat = [dataConcat, data.trial{i}];
    if i == 1
        nextTime = data.time{i};
    elseif i > 1
        nextTime = data.time{i} + max(timeConcat);
    end
    timeConcat = [timeConcat, nextTime];
end

dataConcat(isnan(dataConcat)) = 0; % turn nans into zero

concatData = [];

% concatData.hdr = data.hdr;
concatData.label = data.label;
concatData.fsample = data.fsample;
concatData.trial = {dataConcat};
concatData.time = {timeConcat};
%     concatData.dimord = 'chan_time';
%     concatData.sampleinfo = [1 104381];
%     concatData.elec = ft_read_sens(fileID);
%     concatData.trialinfo = clean_data.trialinfo;
%     concatData.sampleinfo = clean_data.sampleinfo;
%     concatData.cfg = clean_data.cfg;

end