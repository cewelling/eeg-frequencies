%% Analyse Binoc with FFT
clearvars;
clc;
% load('buttons_freqs.mat');
% if exist('E:\Documents\Recorded Data\EEG Feb 2015', 'dir') % location on desktop
%     file_directory = 'E:\Documents\Recorded Data\EEG Feb 2015';
% elseif exist('D:\Recorded Data', 'dir') % location on laptop
%     file_directory = 'D:\Recorded Data';
% elseif exist('C:\EEG Data\mit-data', 'dir')
%     file_directory = 'C:\EEG Data\mit-data';
% elseif exist('C:\EEG Data\mit-data', 'dir')
%     file_directory = 'C:\EEG Data\mit-data';
% else
%     error('please provide directory where file is stored');
% end

%% Clear everything and establish where data is
if exist('E:\Documents\Recorded Data\EEG Feb 2015', 'dir') % location on desktop
    file_directory = 'E:\Documents\Recorded Data\EEG Feb 2015';
elseif exist('D:\Recorded Data', 'dir') % location on laptop
    file_directory = 'D:\Recorded Data';
elseif exist('C:\EEG Data\mit-data', 'dir')
    file_directory = 'C:\EEG Data\mit-data';
elseif exist('/Users/jacksonclee/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results', 'dir')
    file_directory = '/Users/jacksonclee/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results';
elseif exist('/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results', 'dir')
    file_directory = '/Users/robertsonce/Dropbox (MIT)/Projects/EEG_Frequencies/runScripts/Results';
else
    error('please provide directory where file is stored');
end

%% File Names

n{1} = 1;
analRunNum = 1;

%% Analysis Variables
% trial_dur = 30; %JL 12; %JF
discard_start = 3; % how much time should be cut off at beginning %0.5
occipitals = [20 26 27 28 29 30 31];
electrodes = [20 26 28 29 30 31]; %[1:32]; %[27, 29, 30,28]; %occipitals [1:19 21:32]

%stimulation_frequencies = [28.33]; %, 21.25];
% intermodulation_frequencies = [21.6, 43.2];
% intermodulation_frequencies = 21.6;

% Dates & Day #: day2 - 3_11; day3 - 3_14; day4 - 3_15; day5 - 3_17; day6 -
% 3_21; day7 - 3_24

fileNames = dir('../../runScripts/Results/stratus03*3_29*');
fileNames = {fileNames.name};
filenames{1}(1).name = char(fileNames(analRunNum));
% legend = importdata('pilotlegends/pilot-day7-caroline.txt');
legend = importdata('pilotlegends/stratus03.txt');

stimulation_frequencies = legend.data(analRunNum,1:2);
trialNum = legend.data(analRunNum,3);
trial_dur = legend.data(analRunNum,4);


for group = 1:1
    for subject = 1:n{group}
    %% Update the progress bar, load psychometric data
    clc; fprintf('Processing: Group %2.0f of 2, Subject %2.0f of %2.0f\n', group, subject, n{group});
    fileID = fullfile(file_directory, filenames{group}(subject).name);


    %% Trial Definition
    cfg_trldef = [];
    cfg_trldef.dataset = fileID;
    cfg_trldef.trialdef.eventtype = 'STATUS';
    cfg_trldef.trialfun = 'ft_trialfun_general';
    cfg_trldef.trialdef.prestim = -discard_start;
    cfg_trldef.trialdef.poststim = trial_dur; %actual length of trial 12 s
%     cfg_trldef.trialdef.eventvalue = 201:(200+trialNum);
    cfg_trldef.trialdef.eventvalue = 205:208;

    try
        cfg_trldef = ft_definetrial(cfg_trldef);
    catch define_trial_error
        cfg_trldef.trialdef.eventvalue
        cfg_trldef.dataset
        rethrow(define_trial_error);
    end

    %     cfg_trldef.trl = remove_overlaps(cfg_trldef); % this script removes overlapping trials in case triggers got confused

    %% Preprocessing

    cfg_preproc = cfg_trldef;
    cfg_preproc.channel = 'all';
    cfg_preproc.continuous = 'yes';
    cfg_preproc.demean    = 'yes';
    cfg_preproc.detrend = 'yes';
    cfg_preproc.reref = 'yes';
    cfg_preproc.refchannel = 'all';

    cfg_preproc.hpfilter = 'yes';
    cfg_preproc.hpfreq = 2;

    all_data = ft_preprocessing(cfg_preproc);
    
    cfg        = [];
    cfg.method = 'channel';
%     cfg = ft_databrowser(cfg, all_data);
    ft_rejectvisual(cfg, all_data)


    %% FFT
    cfg_fft = [];
    cfg_fft.continuous = 'yes';
    cfg_fft.output = 'pow';
    cfg_fft.method = 'mtmfft';
    cfg_fft.foilim = [5, 30]; %[5, 30];
    cfg_fft.tapsmofrq = 0.09;
    cfg_fft.channel = electrodes; % electrodes or occipitals
    cfg_fft.keeptrials = 'no';

    freqs = ft_freqanalysis(cfg_fft, all_data);
    group_freqs{group, subject} = freqs;


    %% Plot the spectrum

    figure;
    for iTrial = 4%1:1;

    %     subplot(4, 4, iTrial);

        plot(freqs.freq, squeeze( freqs.powspctrm(:, :) ));
        axis([5 30 0 0.2]) %0.18
    %     gridxy([5, 12]);

    %     xlim([3, 10]);
    %     ylim([0, 15]);

    %     gridxy(stimulation_frequencies);
    end


    end % End of participant loop
end % End of group loop








