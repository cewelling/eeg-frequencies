function runRLS(stimulation_frequencies,intermodulation_frequencies, analRunNum, filenames, trialNum, trial_dur,parName,legend,maxIndices)

%% Clear everything and establish where data is
% clearvars;
% clc;

maxIndices

%% Analysis Variables
for electrodeSet = 2 %1:2
    
    %% Load default anal params
    defaultAnalParams
    
    if electrodeSet == 1
        analElectrodes = electrodes;
        electrodeNames = 'allElectrodes';
    elseif electrodeSet == 2
        analElectrodes = occipitals;
        electrodeNames = 'occipitals';
    end


    %% Practice Analysis with one subject
    fileID = fullfile(file_directory, filenames{1}(1).name);

    %% load config params
    [cfg_trldef, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,analElectrodes)

    %% preprocess data
    all_data = ft_preprocessing(cfg_preproc);

    %% Percept Definition
    % Find all Percept Starts
    cfg_percdef.dataset = fileID;
    cfg_percdef.trialdef.eventtype = 'STATUS';
    cfg_percdef.trialfun = 'ft_trialfun_general';
    cfg_percdef.trialdef.prestim = 0;
    cfg_percdef.trialdef.poststim = 0;
    cfg_percdef.trialdef.eventvalue = 1:10;
    cfg_percdef = ft_definetrial(cfg_percdef);


%     for iTrial = 1:trialNum
% 
%         trial_start = cfg_trldef.trl(iTrial, 1);
%         trial_end = cfg_trldef.trl(iTrial, 2);
% 
%         currPercepts = (cfg_percdef.trl(:, 1) >= trial_start) & (cfg_percdef.trl(:, 1) <= trial_end);
% 
%         percepts(iTrial).trl = cfg_percdef.trl(currPercepts, :);
%         percepts(iTrial).type = dec2bin(cfg_percdef.trl(currPercepts, 4)-1) - '0';
%         percepts(iTrial).start = cfg_percdef.trl(currPercepts, 1);
%         percepts(iTrial).duration = [percepts(iTrial).start(2:end); trial_end] - percepts(iTrial).start;
%         percepts(iTrial).start = (percepts(iTrial).start-trial_start) /all_data.fsample;
%         percepts(iTrial).duration = percepts(iTrial).duration /all_data.fsample;
% 
%     end

%% Set up

    % Find and remove artifacts
%     clean_data = removeArtifacts(cfg_trldef,'yes');

%     % Concatenate trials
%     concatData = concatenateTrials(all_data);
%     
%     % Create data files corresponding to perceptual reports
%     [lowFreqReportData, highFreqReportData] = epochOfInt(omniData, concatData);
%     
%     all_data = lowFreqReportData;
%     trial_dur = ceil(max(all_data.time{1,1}));
%     trialNum = 1;

    %% RLS Analysis
    rls_data = all_data;
    single_rls(1) = all_data;
    single_rls(2) = all_data;
    msglength = 0;
    for iTrial = 1:size(all_data.trial, 2);

        fprintf(repmat('\b',1,msglength));
        msg = sprintf('Running RLS - Processing Trial No %d\n', iTrial);
        fprintf(msg);
        msglength = numel(msg);

        colours = {[83 148 255]/255, [255 117 117]/255};
        rls_data.trial{iTrial} = zeros(size(rls_data.trial{iTrial}));

        % Fundamental Stimulation Frequency
        for iFreq = stimulation_frequencies
            j = find(iFreq==stimulation_frequencies);
            cfg_rls.n_cycles = (trial_dur - discard_start) * iFreq;
            cfg_rls.stim_freq = iFreq;
            cfg_rls.channel = focus_electrode;

            [single_rls(j).trial{iTrial}, single_rls(j).amp{iTrial}] = rls_slave( cfg_rls, all_data.trial{iTrial} );

            rls_data.trial{iTrial} = rls_data.trial{iTrial} + single_rls(j).trial{iTrial};
        end


    end

    %% Plot the RLS Signal Amplitude Estimate

    timeFig = figure;
    meanFig = figure;
    % scatterFig = figure;
    
    trialdifferences = []; %store the difference scores
    
    for iElectrode = maxIndices
        for iTrial = 1:size(all_data.trial, 2);

            figure(timeFig);
            subplot(round(trialNum/2), 3, iTrial); hold on;

            for iFreq = 1:2;

                y{iFreq} = smooth(single_rls(iFreq).amp{iTrial}(focus_electrode, :), 150);        
                plot(y{iFreq}, 'Color', colours{iFreq});

            end

            trialdifferences(iTrial,:) = y{1} - y{2};

    %         figure(scatterFig);
    %         
    %         subplot(4, 4, iTrial); hold on; xlim([-4, 4]); ylim([-4, 4]);
    %         time_index = single_rls(1).time{1} > 2;
    %         xscatt = zscore(single_rls(1).amp{iTrial}(focus_electrode, time_index));
    %         yscatt = zscore(single_rls(2).amp{iTrial}(focus_electrode, time_index));
    %         scatter(xscatt, yscatt, 1);
        end
    end



    %% FFT Analysis
    % compare the FFT amplitude between RLS estimate and raw data 
    cfg_fft = [];
    cfg_fft.continuous = 'yes';
    cfg_fft.output = 'pow';
    cfg_fft.method = 'mtmfft';
    cfg_fft.keeptrials = 'yes'; %yes
    cfg_fft.foilim = [0 60]; %[0 60]
    cfg_fft.tapsmofrq = 0.09;
    cfg_fft.channel = 1:32;

    freqs_raw = ft_freqanalysis(cfg_fft, all_data);
    freqs_rls = ft_freqanalysis(cfg_fft, rls_data);

    figure;
    electrode_to_plot = focus_electrode; %focus_electrode
    for iTrial = 1:trialNum
        subplot(round(trialNum/4), 3, iTrial); hold on; % changed round to ceil: round(trialNum/4), 5

        plot(freqs_raw.freq, squeeze( freqs_raw.powspctrm(iTrial, electrode_to_plot, :) ), 'Color', colours{1},'LineWidth',2);

        plot(freqs_rls.freq, squeeze( freqs_rls.powspctrm(iTrial, electrode_to_plot, :) ), '--r');
        xlim([min(stimulation_frequencies)-5, max(stimulation_frequencies)+5]);
    %     ylim([0, 1]);
    %     xlim([0, max([10, stimulation_frequencies])+10]);
    end

    box off
    title(['RLS Fit Results: ' parName ', run ' num2str(analRunNum) ' , ' electrodeNames ', ' num2str(trialNum) ' trials'])
    figName = [plotDir parName '_run' num2str(analRunNum) '_RLSfits_' electrodeNames '.jpg']
    saveas(gcf,figName,'jpg')

    
    %% SNR

    for iFreq = 1:length(stimulation_frequencies)

        snr(iFreq, 1:32) = ...
            freqs_raw.powspctrm(iTrial, :, find(freqs_raw.freq>stimulation_frequencies(iFreq), 1))...
            ./median(freqs_raw.powspctrm(iTrial, :, (freqs_raw.freq>stimulation_frequencies(iFreq)-3 & freqs_raw.freq<stimulation_frequencies(iFreq)+3)), 3);

    end


    % return;



    %% TFR Analysis
    % Trial Definition
    cfg_tfr = [];
    cfg_preproc = [];
    discard_start = 0;
    [cfg_tfr, cfg_preproc] = configParams(fileID,discard_start,trial_dur,trialNum,analElectrodes);
    all_data = ft_preprocessing(cfg_preproc);

    cfg_tfr = [];
    cfg_tfr.channel = analElectrodes; %'all' or analElectrodes
    cfg_tfr.foi = 0:(1/30):maxFreq; %frequency of interest  %original: step by 0.25
    cfg_tfr.toi = -3:slidingWindowStep:trial_dur; %time of interest
    cfg_tfr.keeptrials = 'yes'; % this makes sure all trials are separate

    cfg_tfr.method = 'mtmconvol'; % cfg_tfr.method = 'wavelet';
    cfg_tfr.taper = 'hanning'; %sliding window shape (tapered)
    cfg_tfr.t_ftimwin = slidingWindowWidth*ones(size(cfg_tfr.foi)); % cfg_tfr.width = 16;

    tfr_raw = ft_freqanalysis(cfg_tfr, all_data);
    tfr_rls = ft_freqanalysis(cfg_tfr, rls_data);
    
    %% Save RLS data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataNameRAW = [rlsDir parName '_RLSRAW_run' num2str(analRunNum) '_' electrodeNames '.mat'];
    dataNameRLS = [rlsDir parName '_RLSModel_run' num2str(analRunNum) '_' electrodeNames '.mat'];
    
    save(dataNameRAW,'tfr_raw')
    save(dataNameRLS,'tfr_rls')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    % for iFreq = 1:2
    %  tfr_single_freq{iFreq} = ft_freqanalysis(cfg_tfr, single_rls(iFreq));
    % end
    cmap = cbrewer('seq', 'PuBu', 200);

    yLims = [1 0.8; 1 1; 1 1]; %trial 3
    zLims = [0.01 .1; 0.01 .07; 0.01 0.07]; %trial 5
    for iTrial = 1:trialNum
        figure;
        for iFreq = 1:length(stimulation_frequencies)
            subplot(2, 1, iFreq);
            cfg_plot = [];
            cfg_plot.channel = freqs_raw.label(focus_electrode);
            cfg_plot.channel = 'all';
            cfg_plot.interactive = 'no';
            cfg_plot.trials = iTrial;
            cfg_plot.xlim = [0 trial_dur];
            cfg_plot.ylim = [stimulation_frequencies(iFreq)-yLims(iFreq,1), stimulation_frequencies(iFreq)+yLims(iFreq,2)];
            %cfg_plot.zlim = [zLims(iFreq,1),zLims(iFreq,2)];
            %cfg_plot.baseline = [-3, -1];
            cfg_plot.baselinetype = 'relative';

            gridxy2([], stimulation_frequencies);

            ft_singleplotTFR(cfg_plot, tfr_rls  );
            %colormap(cmap);
        end
        box off
        title(['RLS Heat Results: ' parName ', run ' num2str(analRunNum) ' , ' electrodeNames ', ' num2str(trialNum) ' trials'])
        figName = [plotDir parName '_run' num2str(analRunNum) '_RLSheat_' electrodeNames '.jpg'];
        saveas(gcf,figName,'jpg')

    end

  clearvars -except *frequencies  analRunNum  filenames  trialNum  trial_dur parName legend
end
