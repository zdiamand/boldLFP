
% addpath 'C:\Users\18133\OneDrive\Documents\MATLAB\fieldtrip-20220707'
% ft_defaults

load('D:\Olab\patientData\extractedNoCut\P62CS_041919_lfpLoopNoCut.mat');

% Define a cell array to store the data for each channel
num_channels = length(dat);
time_series_data_channels = cell(1, num_channels);

% Store data from each channel in a new array
for i = 1:num_channels
    time_series_data_channels{i} = dat{i}.dataRawDownsampled;  %1x16 cell of 1x1 cells of Nx1 doubles
end

% Concatenate all of the channels into a single array
time_series_data = cat(2, time_series_data_channels{:}); %1x16 cell of Nx1 doubles
time_series_data_matrix = transpose(cell2mat(time_series_data)); %16xN cell of 
time_series_data_matrix2 = {time_series_data_matrix};

% Timepoints should be the same for all channels, so we use the first
timepoints = dat{1}.timestampDownsampled;

% Convert to seconds 
timepoints = cellfun(@(d) arrayfun(@(x) x / 10^6, d), timepoints, 'UniformOutput', false);

% Trial lengths for later use
trial_lengths = (dat{1}.trialEndTTLs - dat{1}.trialStartTTLs) /10^6;


channelLabels = cell(num_channels, 1); % initialize the cell array
for i = 1:num_channels
    channelLabels{i, 1} = ['channel' num2str(i)]; % store the strings in the cell array
end


% Data setup in proper format
data = [];
data.label = channelLabels; % channel label
data.time = timepoints; 
data.trial = time_series_data_matrix2; % time series data
data.fsample = 1000; % sampling frequency


cfg = [];
cfg.artfctdef.zvalue.channel = channelLabels;
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.cutoff = 20;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.cumulative = 'yes';
cfg.artfctdef.zvalue.medianfilter = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff = 'yes';

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';


[cfg, artifact_jump] = ft_artifact_zvalue(cfg, data);


cfg                           = [];
cfg.artfctdef.reject          = 'partial'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
cfg.artfctdef.jump.artifact   = artifact_jump;
data_no_artifacts = ft_rejectartifact(cfg,data);



% Perform ICA
num_components = 10;

cfg = [];
cfg.method = 'runica'; % ICA method to use, in this case runica
%running with 10 components, because this subject only has 16 channels
cfg.numcomponent = num_components; % number of components to retain
comp = ft_componentanalysis(cfg, data_no_artifacts);

kl_divergences = zeros(1, num_components); % initialize an array to store the KL divergences
for i = 1:num_components
    columndata = comp.topo(:, i); % extract the data from the i-th column
    pmf_data = histcounts(columndata, num_channels, 'Normalization', 'probability'); % estimate the PMF of the data
    pmf_data(pmf_data == 0) = 0.00001; % add a small positive value to any element equal to 0
    pmf_uniform = ones(1, num_channels) / num_channels; % PMF of the uniform distribution
    kl_divergences(i) = sum(pmf_data .* log(pmf_data ./ pmf_uniform)); % calculate the KL divergence
end

figure;
bar(kl_divergences);
title('KL Divergences for Each Component');
xlabel('Component');
ylabel('KL Divergence');


threshold = 0.8; % The threshold for what is considered "equal distribution"
equal_components = find(kl_divergences < threshold);

% Display the indices of the equal components
cfg = [];
cfg.component = equal_components; % to be removed
data_clean = ft_rejectcomponent(cfg, comp);


% Perform time-frequency analysis using the mtmconvol method
cfg = [];
cfg.output = 'pow';
cfg.channel = channelLabels;
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 2:1:30; %frequency of interest
cfg.t_ftimwin = 7./cfg.foi;  % 7 cycles per time window %NEED TO CHANGE
cfg.toi = min(cell2mat(timepoints)):0.05:max(cell2mat(timepoints));
tf_mtmconvol = ft_freqanalysis(cfg, data_clean);



% Load sliced version of data to extract the trial start and end times for
% each trial 
load('D:\Olab\patientData\extracted\P62CS_041919_lfpLoop.mat');
trialStartTTLs = dat{1}.trialStartTTLs / 10^6;
trialEndTTLs = dat{1}.trialEndTTLs / 10^6;
trialLengths = trialEndTTLs - trialStartTTLs;
numTrials = length(trialStartTTLs);
trialData = cell(1, numTrials);
time_trials = cell(numTrials, 1);

%Cut frequency data into trial slices
for i = 1:numTrials
    startIndex = find(tf_mtmconvol.time >= trialStartTTLs(i), 1);
    endIndex = find(tf_mtmconvol.time <= trialEndTTLs(i), 1, 'last');
    trialData{i} = tf_mtmconvol.powspctrm(:, :, startIndex:endIndex);
    time_trials{i} = tf_mtmconvol.time(startIndex:endIndex);
end

% Shift the timepoints so they are from -1 to end of trial
result_cell = cell(1, length(time_trials));
    for i = 1:length(time_trials)
    d = time_trials{i};
    result_cell{i} = shift_timepoints(d);
    end
timepoints = result_cell;


% Put each trial's data into correct data structure
nTrials = numel(trialData);
freq = cell(1, nTrials);
for i = 1:nTrials
    freq{i}.powspctrm = cell2mat(trialData(i));
    freq{i}.label = tf_mtmconvol.label;
    freq{i}.freq = tf_mtmconvol.freq;
    freq{i}.time = cell2mat(timepoints(i));
    freq{i}.dimord = 'chan_freq_time';
end


%Average over all trials. This will cut all trials to the shortest trial's
%length
cfg = [];
cfg.keepindividual = 'no';
zz = ft_freqgrandaverage(cfg, freq{:});

%Plot with no baseline correction
ft_singleplotTFR(cfg, zz);

    
% Plot the results with baseline correction
cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'relative';
cfg.maskstyle    = 'saturation';
cfg.channel = channelLabels;
cfg.interactive  = 'no';
figure
ft_singleplotTFR(cfg, zz);
