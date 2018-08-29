
for electrodeSet = 2%1:2
    
    if electrodeSet == 1
        analElectrodes = electrodes;
        electrodeNames = 'allElectrodes';
    else
        analElectrodes = occipitals;
        electrodeNames = 'occipitals';
    end
    dataNameRAW = [rlsDir parName '_RLSRAW_run' num2str(analRunNum) '_' electrodeNames '.mat'];
    dataNameRLS = [rlsDir parName '_RLSModel_run' num2str(analRunNum) '_' electrodeNames '.mat'];
end

if exist(dataNameRLS)
    %load dataNameRAW
    load(dataNameRLS) %perform analysis on RLS data
else
    sprintf('Error: You need to run the RLS analysis first!')
end

if sum(round(stimulation_frequencies) == round([14.16 17])) == 2
    if length(tfr_rls.freq) == 120
        channelIndices1 = 56:58;
        channelIndices2 = 67:69;
        interFreqIndices = 44:46; %11.33 Hz
    elseif length(tfr_rls.freq) == 900
        channelIndices1 = 425;
        channelIndices2 = 510;
        interFreqIndices = 340;
    end
elseif sum(round(stimulation_frequencies) == round([17 21.25])) == 2
    channelIndices1 = 67:69;
    channelIndices2 = 84:86;
    interFreqIndices = 50:52; %12.75 Hz
end

