
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


channelLabels = cell(num_channels, 1); % initialize the cell array
for i = 1:num_channels
    channelLabels{i, 1} = ['channel' num2str(i)]; % store the strings in the cell array
end


% Data setup in proper format
data = [];
data.label = channelLabels; % channel label
data.time = timepoints; 
data.trial = time_series_data; % time series data
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


cfg = [];
cfg.artfctdef.zvalue.channel = channelLabels;
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.cutoff       = 4;
cfg.artfctdef.zvalue.trlpadding   = 0;
cfg.artfctdef.zvalue.fltpadding   = 0;
cfg.artfctdef.zvalue.artpadding   = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter     = 'yes';
cfg.artfctdef.zvalue.bpfreq       = [110 140];
cfg.artfctdef.zvalue.bpfiltord    = 9;
cfg.artfctdef.zvalue.bpfilttype   = 'but';
cfg.artfctdef.zvalue.hilbert      = 'yes';
cfg.artfctdef.zvalue.boxcar       = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_muscle] = ft_artifact_zvalue(cfg, data);

%cfg = [];
%cfg.artfctdef.zvalue.channel = channelLabels;
%cfg.artfctdef.zvalue.cutoff      = 4;
%cfg.artfctdef.zvalue.trlpadding  = 0;
%cfg.artfctdef.zvalue.artpadding  = 0.1;
%cfg.artfctdef.zvalue.fltpadding  = 0;

% algorithmic parameters
%cfg.artfctdef.zvalue.bpfilter   = 'yes';
%cfg.artfctdef.zvalue.bpfilttype = 'but';
%cfg.artfctdef.zvalue.bpfreq     = [2 15];
%cfg.artfctdef.zvalue.bpfiltord  = 4;
%cfg.artfctdef.zvalue.hilbert    = 'yes';

% feedback
%cfg.artfctdef.zvalue.interactive = 'yes';

%[cfg, artifact_eog] = ft_artifact_zvalue(cfg, data);






cfg                           = [];
cfg.artfctdef.reject          = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
cfg.artfctdef.jump.artifact   = artifact_jump;
cfg.artfctdef.muscle.artifact = artifact_muscle;
%cfg.artfctdef.eog.artifact    = artifact_eog;
data_no_artifacts = ft_rejectartifact(cfg,data);



% Perform ICA
num_components = 10;

cfg = [];
cfg.method = 'runica'; % ICA method to use, in this case runica
%running with 10 components, because this subject only has 16 channels
cfg.numcomponent = num_components; % number of components to retain
comp = ft_componentanalysis(cfg, data_no_artifacts);

% Calculate the mean variacne of each components spatial distribution
%var_magnitudes = var(abs(comp.topo));

kl_divergences = zeros(1, num_components); % initialize an array to store the KL divergences
for i = 1:num_components
    columndata = comp.topo(:, i); % extract the data from the i-th column
    pmf_data = histcounts(columndata, num_channels, 'Normalization', 'probability'); % estimate the PMF of the data
    pmf_data(pmf_data == 0) = 0.00001; % add a small positive value to any element equal to 0
    pmf_uniform = ones(1, num_channels) / num_channels; % PMF of the uniform distribution
    kl_divergences(i) = sum(pmf_data .* log(pmf_data ./ pmf_uniform)); % calculate the KL divergence
end


threshold = 1; % The threshold for what is considered "equal distribution"
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
cfg.t_ftimwin = 7./cfg.foi;  % 7 cycles per time window
cfg.toi = -1:0.05:max(trial_lengths); %from -1 seconds to the max trial length
tf_mtmconvol = ft_freqanalysis(cfg, data_clean);
    
% Plot the results
cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'relative';
cfg.maskstyle    = 'saturation';
cfg.channel = channelLabels;
cfg.interactive  = 'no';
figure
ft_singleplotTFR(cfg, tf_mtmconvol);
