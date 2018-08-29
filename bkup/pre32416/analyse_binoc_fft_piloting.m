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

file_directory = '/Users/jacksonclee/Documents/MATLAB/Rivalry/data/';

%% File Names
% filenames{1} = dir([file_directory, '\*CTR*BinSSVEP.bdf']);
% filenames{2} = dir([file_directory, '\*ASC*BinSSVEP.bdf']);
% n{1} = size(filenames{1}, 1);
% n{2} = size(filenames{2}, 1);

n{1} = 1;
filenames{1}(1).name = 'jl_2_100216.bdf';

%% Analysis Variables
trial_dur = 12;
discard_start = 0.5; % how much time should be cut off at beginning
% electrodes = [27, 29, 64];
electrodes = [1:32];

stimulation_frequencies = [5, 8.5];
% intermodulation_frequencies = [21.6, 43.2];
% intermodulation_frequencies = 21.6;


% group_freqs = cell(2, max([n{1}, n{2}]));


for group = 1:1
for subject = 1:n{group}
%% Update the progress bar, load psychometric data
    clc; fprintf('Processing: Group %2.0f of 2, Subject %2.0f of %2.0f\n', group, subject, n{group});
    
    fileID = fullfile(file_directory, filenames{group}(subject).name);
%     psychometric_fileID = [fileID(1:end-13), '.mat'];
%     official_ID = str2double(fileID(end-16:end-13));
%     try
%         psycho_data = load(psychometric_fileID);
%         aq{group}(subject, 1) = psycho_data.AQ;
%     catch
%         aq{group}(subject, 1) = NaN;
%     end
    
    
%% Trial Definition
    cfg_trldef = [];
    cfg_trldef.dataset = fileID;
    cfg_trldef.trialdef.eventtype = 'STATUS';
    cfg_trldef.trialfun = 'ft_trialfun_general';
    cfg_trldef.trialdef.prestim = -discard_start;
    cfg_trldef.trialdef.poststim = trial_dur; %actual length of trial 12 s
    cfg_trldef.trialdef.eventvalue = 201:216;

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
cfg_preproc.channel = 1:32;
cfg_preproc.continuous = 'yes';
cfg_preproc.demean    = 'yes';
cfg_preproc.detrend = 'no';
cfg_preproc.reref = 'yes';
cfg_preproc.refchannel = electrodes;

cfg_preproc.hpfilter = 'yes';
cfg_preproc.hpfreq = 2;

% cfg_preproc.refchannel = [27 28 64 30];

all_data = ft_preprocessing(cfg_preproc);

%% FFT
cfg_fft = [];
cfg_fft.continuous = 'yes';
cfg_fft.output = 'pow';
cfg_fft.method = 'mtmfft';
cfg_fft.foilim = [5, 45];
cfg_fft.tapsmofrq = 0.09;
cfg_fft.channel = electrodes;
cfg_fft.keeptrials = 'no';

freqs = ft_freqanalysis(cfg_fft, all_data);
group_freqs{group, subject} = freqs;


%% Plot the spectrum

figure;
for iTrial = 1:1;
    
%     subplot(4, 4, iTrial);
    
    plot(freqs.freq, squeeze( freqs.powspctrm(:, :) ));
    
%     gridxy([5, 12]);
    
%     xlim([3, 10]);
%     ylim([0, 15]);
    
%     gridxy(stimulation_frequencies);
end


end % End of participant loop
end % End of group loop








