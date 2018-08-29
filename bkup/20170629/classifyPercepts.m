% Uses classifier to...
% uses the Neural Decoding Toolbox (http://www.readout.info/)

clearvars

% load parameters
analysisParams

% participant to test
parName = 'stratus137';
latency = round(.8*sampRate); % predicted latency in data points
segPoints = round(0.5*sampRate); % segment length in data points (each segment is a "trial" to be classified)

%% Put data in raster format

% create and save rasters
makeRaster

%% bin the data 
bin_width = 1; % # time points
step_size = 1; % # time points
saveName = 'ML/binned_data';
create_binned_data_from_raster_data('ML/raster_data/', saveName, bin_width, step_size); 

%% create a datasource (DS) object
binnedFile = 'ML/binned_data_1ms_bins_1ms_sampled.mat';
labelToDecode = 'domFreq';
numSplits = 2; % for numSplits = 5, 4/5 of data become training set

% differentiate training from testing set
trainLabels{1} = {'lo'};
testLabels{1} = {'lo_t'};
trainLabels{2} = {'hi'};
testLabels{2} = {'hi_t'};

ds = generalization_DS(binnedFile, labelToDecode, numSplits, trainLabels, testLabels);

% Find number of percept repeats
load(['ML/raster_data/' parName '_lowFreq.mat']);
hiRepeats = length(find(strcmp(raster_labels.domFreq,'hi')));
loRepeats = length(find(strcmp(raster_labels.domFreq,'lo')));

% Use max number of percept repeats possible per testing / training set
ds.num_times_to_repeat_each_label_per_cv_split = min([hiRepeats loRepeats]) / numSplits;

%%
% create a preprocessor
% zscore_normalize? Otherwise higher frequency will contribute more to the
% model than lower frequency
%preprocessor = ;

%% create a classifier
classifier = max_correlation_coefficient_CL;

%% create a cross-validator
xValidator = standard_resample_CV(ds, classifier); 
xValidator.num_resample_runs = 1;

%% run the decoding analysis
DECODING_RESULTS = xValidator.run_cv_decoding;
%%
saveName = ['ML/' parName '_simpleTest'];
save(saveName, 'DECODING_RESULTS');

%% plot decoding accuracy as a function of time
resultNames{1} = saveName;
plot_obj = plot_standard_results_object(resultNames);
plot_obj.significant_event_times = 0; % start of reported epoch
plot_obj.plot_results;