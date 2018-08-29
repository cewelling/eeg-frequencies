%%fftParams: called by runFFT and runRLS
%%Contains all params for preprocessing the raw data and running the FFT
%%Note: re-referencing and bandpass params are in ftPreProc

if strcmp(parName, 'cumulus25')
    elecsUsed = 33:64;
else
    elecsUsed = 1:32;
end

%Common preprocessing params
cfg_preproc.channel = elecsUsed;
cfg_preproc.continuous = 'yes';
cfg_preproc.demean    = 'no'; % perform baseline correction
cfg_preproc.detrend = 'yes'; % remove linear trend from the data
cfg_preproc.bpfilter = 'yes'; % bandpass filter?
cfg_preproc.bsfilter = 'yes'; % bandstop filter?
cfg_preproc.bsfreq = [59 61]; % bandstop frequencies (US = 60)

%Params for running FFT after preprocessing
analElectrodes = electrodeSet.nums;
cfg_fft = [];
cfg_fft.output = 'pow'; % returns the power-spectra
cfg_fft.method = 'mtmfft'; % analyzes entire spectrum
cfg_fft.foilim = [5, 30]; %[5, 30];
cfg_fft.taper = 'hanning';
% cfg_fft.tapsmofrq = .05; % multitaper spectral smoothing in Hz
cfg_fft.channel = analElectrodes; % electrodes or occipitals
cfg_fft.keeptrials = 'yes'; % return individual trials, not average
% cfg_fft.pad = ceil(max(concatData.time{1,1})); % might be useful
% if we want to compare spectra between groups
