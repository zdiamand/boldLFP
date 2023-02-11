load('D:\Olab\patientData\extracted\P62CS_041919_lfpLoop.mat');

% addpath 'C:\Users\18133\OneDrive\Documents\MATLAB\fieldtrip-20220707'
% ft_defaults

% Note: only using two channels for now to get the code right

% Each cell is a d_i x 1 double. Need to zero pad so that we have d x 1
time_series_data_channel1 = pad_cell_array(dat{1}.dataRawDownsampled);
time_series_data_channel2 = pad_cell_array(dat{2}.dataRawDownsampled);

% Join together so each cell  (1 x num trials) is in shape (num channels x num samples) 
time_series_data = join_cell_arrays(time_series_data_channel1, time_series_data_channel2);

% Each cell is a 1 x d_i double. Need to zero pad and transpose so that we have d x 1
% Timepoints should be the same for all channels, so we use the first
timepoints = pad_cell_array(dat{1}.timestampDownsampled);

% Convert to seconds 
timepoints = cellfun(@(d) arrayfun(@(x) x / 10^6, d), timepoints, 'UniformOutput', false);

% Trial lengths for later use
trial_lengths = (dat{1}.trialEndTTLs - dat{1}.trialStartTTLs) /10^6;

% Shift the timepoints so they are from -1 to end 
result_cell = cell(1, 300);
    for i = 1:300
    d = timepoints{i};
    result_cell{i} = shift_timepoints(d);
    end
timepoints = result_cell;

% Data setup in proper format
data = [];
data.label = {'channel1'; 'channel2'}; % channel label
data.time = timepoints; 
data.trial = time_series_data; % time series data
data.fsample = 1000; % sampling frequency

% Preprocess the data (e.g. filtering, resampling, etc.)
% TODO: Artifact removal 
data_preprocessed = data; 
%data_preprocessed = ft_preprocessing(cfg, data);

% Perform time-frequency analysis using the mtmconvol method
cfg = [];
cfg.output = 'pow';
cfg.channel = {'channel1'; 'channel2'};
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 2:1:30; %frequency of interest
cfg.t_ftimwin = 7./cfg.foi;  % 7 cycles per time window
cfg.toi = -1:0.05:max(trial_lengths); %from -1 seconds to the max trial length
tf_mtmconvol = ft_freqanalysis(cfg, data_preprocessed);
    
% Plot the results
cfg              = [];
cfg.baseline     = [-1 max(trial_lengths)];
cfg.baselinetype = 'relative';
cfg.maskstyle    = 'saturation';
%cfg.zlim         = [-2e-27 2e-27];
cfg.channel = {'channel1'; 'channel2'};
cfg.interactive  = 'no';
figure
ft_singleplotTFR(cfg, tf_mtmconvol);
