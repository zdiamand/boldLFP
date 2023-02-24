% This function loops macro channels and extracts the raw LFP signal
function data = loopLFPChannels_standalone(channelIds,session)
dataReplay = fullfile('D:\Olab\patientData\', session, '\dataReplay');
contents = dir(dataReplay);
% Get the first subdirectory (ignoring '.' and '..')
for i = 1:length(contents)
    if contents(i).isdir && ~strcmp(contents(i).name, '.') && ~strcmp(contents(i).name, '..')
        rawFolder = fullfile(dataReplay, contents(i).name);
        rawFolder = [rawFolder '\'];
        break;
    end
end

%'D:\Olab\patientData\' patientData '\dataReplay';
sortFolder = [dataReplay '\sort\final\']; %ignore
data = {};
% Extraction mode for CSC files
% extractionMode = 2;

% Apply 60 Hz notch filter?
applyNotchFilter = 1;
% Frequency to downsample 
FsDownReq = 1000;
% Replace spike times with cubic interpolation
replaceMethod = 2; % cubic method
% Time padding before and after trial
padding = 1e6; % 1s
% Looping over available macro channels
for cI = channelIds
    filename = [rawFolder 'CSC' num2str(cI) '.ncs'];
    cell_filename = [sortFolder 'A' num2str(cI) '_sorted_new.mat'];
    allSpikesToRemove = [];
    if exist(cell_filename, 'file') == 2
        spike_data = load(cell_filename);
        allSpikesToRemove = spike_data.newTimestampsNegative;
    end
    if exist(filename, 'file') == 2
        info = struct();        
        [timestamps,nrBlocks,nrSamples,sampleFreq,isContinuous,headerInfo] = ...
            getRawCSCTimestamps(filename);        
        % Getting data from each trial separately
        path = ['D:\Olab\shared\shared\neural\' session '\sessionData.mat'];
        sessiondat = load(path);
        trialStartTTLs = sessiondat.('sessionData').('trialStartTime');
        trialLen = sessiondat.('sessionData').('trialEndTime') * 10^6;
        trialEndTTLs = trialStartTTLs + trialLen;
        nTrials = length(trialStartTTLs);

        trialStartTTLs = min(trialStartTTLs);
        trialEndTTLs = max(trialEndTTLs);

        trialStartTTLs = reshape(trialStartTTLs, [length(trialStartTTLs), 1]);
        trialEndTTLs = reshape(trialEndTTLs, [length(trialEndTTLs), 1]);


        % Extracting portion pertaining to this trial
        [~, dataRawDownsampled, timestampDownsampled, ...
            isContinuous,periodsExtracted, ...
            validTrials, FsDown, Fs, maxRange] = ... 
            getLFPofTrial_casino(filename, trialStartTTLs(:,1), trialEndTTLs(:,1), ...
            applyNotchFilter, FsDownReq, allSpikesToRemove,replaceMethod);
    
        %% Saving channel info to cell
        info.filename = filename;
        info.timestamps = timestamps;
        info.nrBlocks = nrBlocks;
        info.nrSamples = nrSamples;
        info.sampleFreq = sampleFreq;
        info.isContinuous = isContinuous;
        info.headerInfo = headerInfo;
        info.trialStartTTLs = trialStartTTLs;
        info.trialEndTTLs = trialEndTTLs;
        info.nTrials = nTrials;
        info.dataRawDownsampled = dataRawDownsampled;
        info.timestampDownsampled = timestampDownsampled;
        info.isContinuous = isContinuous;
        info.periodsExtracted = periodsExtracted;        
        info.validTrials = validTrials;
        info.FsDown = FsDown;
        info.Fs = Fs;
        info.maxRange = maxRange;
        info.channel = cI;
        info.brainArea = getBrainArea(session,cI);
        data = [data; info];
    end    
    
end
end