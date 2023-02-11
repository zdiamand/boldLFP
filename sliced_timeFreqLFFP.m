load('D:\Olab\patientData\extracted\P62CS_041919_lfpLoop.mat');

% addpath 'C:\Users\18133\OneDrive\Documents\MATLAB\fieldtrip-20220707'
% ft_defaults

% Define a cell array to store the zero-padded data for each channel
num_channels = length(dat);
time_series_data_channels = cell(1, num_channels);

% Zero-pad each channels data and store it in the time_series_data_channels cell array
for i = 1:num_channels
    time_series_data_channels{i} = pad_cell_array(dat{i}.dataRawDownsampled);
end

% Concatenate all of the zero-padded channels into a single array
time_series_data = join_cell_arrays(time_series_data_channels{:});

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
data.label = {'channel1'; 'channel2'; 'channel3'; 'channel4'; 'channel5'; 'channel6'; 'channel7'; 'channel8'; 'channel9'; 'channel10'; 'channel11'; 'channel12'; 'channel13'; 'channel14'; 'channel15'; 'channel16'}; % channel label
data.time = timepoints; 
data.trial = time_series_data; % time series data
data.fsample = 1000; % sampling frequency

% Preprocess the data (e.g. filtering, resampling, etc.)

% TODO: Artifact removal 
% Perform ICA
cfg = [];
cfg.method = 'runica'; % ICA method to use, in this case runica
%running with 10 components, because this subject only has 16 channels
cfg.numcomponent = 10; % number of components to retain
comp = ft_componentanalysis(cfg, data);

% Calculate the mean variacne of each components spatial distribution
var_magnitudes = var(abs(comp.topo));


threshold = 0.001; % The threshold for what is considered "equal distribution"
equal_components = find(var_magnitudes < threshold);

% Display the indices of the equal components
cfg = [];
cfg.component = equal_components; % to be removed
data_clean = ft_rejectcomponent(cfg, comp);

% Perform time-frequency analysis using the mtmconvol method
cfg = [];
cfg.output = 'pow';
cfg.channel = {'channel1'; 'channel2'};
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 2:1:30; %frequency of interest
cfg.t_ftimwin = 7./cfg.foi;  % 7 cycles per time window
cfg.toi = -1:0.05:max(trial_lengths); %from -1 seconds to the max trial length
tf_mtmconvol = ft_freqanalysis(cfg, data_clean);
    
% Plot the results
cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'relative';
cfg.maskstyle    = 'saturation';
cfg.channel = {'channel1'; 'channel2'};
cfg.interactive  = 'no';
figure
ft_singleplotTFR(cfg, tf_mtmconvol);
