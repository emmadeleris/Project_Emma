% We are only working on the no TMR part of the data!

%% Set paths and code variables
clear all

base_dir = "/Users/emma/Desktop/München/MSc Learning Sciences/RP Schreiner/Data analysis";
subjectsFolder = fullfile(base_dir, "/subjects/s*");
fieldtrip_path = "toolbox/fieldtrip-20240916";
breathmetrics_path = "toolbox/breathmetrics";
circstat_path = "toolbox/CircStat2012a";
subjects = dir(subjectsFolder);

fs = 512;
new_fs = 256;

toolbox_paths = {fieldtrip_path, breathmetrics_path, circstat_path};

addpath("/Users/emma/Downloads/FDR.m");

for i = 1:length(toolbox_paths)
    toolbox_path = toolbox_paths{i};

    toolbox_path = fullfile(base_dir, toolbox_path);
    addpath(toolbox_path);

    files = dir(toolbox_path);
    subFolders = {files([files.isdir]).name};
    subFolders = subFolders(~ismember(subFolders, {'.', '..'}));
    fullSubFolders = fullfile(toolbox_path, subFolders);
    for j = 1:length(fullSubFolders)
        addpath(fullSubFolders{j});  % Add each subfolder to the path
    end

end

% Create an empty table with predefined column names and types
subjectsTable = table();

% Predefine the variable names for subjectsTable
eventTypes = {'SO', 'Spindle', 'SO_Spindle'};
channels = {'F1', 'F2', 'F3', 'F5', 'Fz'};
statistics = {'Mean', 'pVal', 'zVal'};

% Create a cell array to hold column names
columnNames = {};
for ievent = 1:length(eventTypes)
    for ichann = 1:length(channels)
        for istat = 1:length(statistics)
            columnNames{end+1} = [eventTypes{ievent}, '_', channels{ichann}, '_', statistics{istat}];
        end
    end
end

% Add mean phase and mean p-value columns for each event type
for ievent = 1:length(eventTypes)
    columnNames{end+1} = [eventTypes{ievent}, '_Mean_Phase'];
    columnNames{end+1} = [eventTypes{ievent}, '_Mean_pVal'];
    columnNames{end+1} = [eventTypes{ievent}, '_Mean_pVal_FDR']; % FDR-corrected p-value
end

% Define StageDurations as a nested table (stored as a cell inside the main table)
StageDurationsTemplate = array2table(zeros(1,5), 'VariableNames', {'Wake', 'N1', 'N2', 'N3', 'REM'});

% Add subject-specific metadata and extra columns
columnNames = [{'SubjectID', 'Age', 'Gender'}, columnNames, {'TST'}, StageDurationsTemplate.Properties.VariableNames];
% Esteban will add the Age and Gender data himself as I don't have access
% to it

% Initialise the main subjectsTable
subjectsTable = cell2table(cell(0, length(columnNames)), 'VariableNames', columnNames);

%% For loop

for isubject = 1:height(subjects)
    isubjectFolder = fullfile(subjects(isubject).folder, subjects(isubject).name);

    signalFile.folder = fullfile(isubjectFolder, "night");
    signalFile.name = 'night_noTMR_Segment_0.edf';

    % Extract Subject ID (last part of folder path)
    subjectID = subjects(isubject).name;

    %% Load the data
    cfg = [];
    cfg.dataset = fullfile(signalFile.folder, signalFile.name);
    data = ft_preprocessing(cfg);

    data.label = cellfun(@(x) x(5:end-4), data.label, 'UniformOutput', false); % Fix label naming

    %% Preprocessing (e.g. changing reference to M1 and M2)
    % VERY IMPORTANT: Preprocess EEG and Respiration separately

    % For EEG:
    cfg             = [];
    cfg.channel         = {'M1' 'M2' 'F1' 'F2' 'F3' 'F5' 'Fz'};
    cfg.reref           = 'yes';
    cfg.refchannel      = {'M1', 'M2'};
    dataEEG = ft_preprocessing(cfg, data);

    % For respiration:
    % Select respiration signal
    cfg = [];
    cfg.channel = 'BIP7'; % check the name!
    dataResp = ft_selectdata(cfg, data);

    dataResp.label{1,1} = 'Resp'; % Rename channel BIP7 to Resp
    dataResp.trial{1,1} = dataResp.trial{1,1}*1e6; % Scale the signal from Volts to uV

    % Append EEG and Resp data
    cfg = [];
    data_all = ft_appenddata(cfg, dataEEG, dataResp);

    % Calculate duration of sleep stages
    hypnoFile = dir(fullfile(signalFile.folder, "hypno_noTMR.txt"));
    hypnoData = load(fullfile(hypnoFile.folder, hypnoFile.name));
    sleepStages = hypnoData(:,1);

    stageCodes = [0, 1, 2, 3, 5]; % Wake, N1, N2, N3, REM
    stageDurations = array2table(arrayfun(@(s) sum(sleepStages == s) * 30 / 60, stageCodes), ...
        'VariableNames', {'Wake', 'N1', 'N2', 'N3', 'REM'});

    % Downsampling
    cfg             = [];
    cfg.resamplefs  = new_fs;

    data_all = ft_resampledata(cfg, data_all);
    data_all.sampleinfo = size(data_all.time{1,1}); % Update sampleinfo to reflect the number of time points in the first trial after resampling

    % Load hypnogram
    cfg              = [];
    cfg.hypnogram    = fullfile(hypnoFile.folder, hypnoFile.name);
    cfg.type         = 'plain'; % import hypno without scored arousal
    cfg.slStLen      = 30; % set the sleep stage length (or epoch duration) to 30 seconds – to properly align the hypnogramme data with the electrophysiological data in data_all
    data_all         = slythm.importHypnogram(cfg, data_all); % create a folder (+slythm) and put the importHypnogram function in it and import it

    % This part is added newly to downsample dataResp to match the dataResp
    % shape with data_all and dataEEG
    cfg = [];
    cfg.channel = {'Resp'};
    dataResp = ft_selectdata(cfg, data_all);

    %% Detect events IN THE EEG (specific channels) - check detecting algos from Esteban

    % For EEG:
    cfg             = [];
    cfg.channel         = {'F1' 'F2' 'F3' 'F5' 'Fz'};
    dataEEG = ft_selectdata(cfg, data_all);

    % Perform SO detection
    cfg                     = [];
    cfg.bpfreq              = [0.3 1.25]; % Sets a band-pass filter with a frequency range of 0.3 Hz to 1.25 Hz, which is the typical frequency range for SOs
    cfg.slStages            = [2 3]; % Specifies that the detection is limited to sleep stages 2 and 3

    cfg.thresCh             = dataEEG.label; % specifies the EEG channels (all of those listed in dataEEG) that will be used for thresholding in the detection process (the thresholding step typically involves defining a cutoff value for detecting significant events)
    cfg.detectCh            = dataEEG.label; % defines the channels where event detection will occur
    cfg.thresType           = 'channelwise'; % the threshold for detecting SOs will be set individually for each channel, rather than using a global threshold across all channels
    cfg.filtType            = 'but'; % specifies the filter type as a Butterworth filter - commonly used in signal processing due to its flat frequency response in the passband

    cfg.criterionLen        = [0.8, 2]; % sets the criteria for the duration of an SO, in seconds
    cfg.criterionParam      = 'both'; % detection criteria is based on both the duration and amplitude of the SO
    cfg.criterionCenter     = 'mean'; % the threshold is based on the mean amplitude
    cfg.criterionVar        = 'scaleCenter'; % scaling factor is applied to determine variability for detection - we adjust the detection threshold to be 1.5 times the standard deviation (or other variability measures) of the signal, effectively making the detection criteria more sensitive to variations in the signal
    cfg.criterionScale      = 1.5; % scaling factor is 1.5, which is often considered a good balance for many types of biological signals
    cfg.criterionPadding    = [0, 0]; % padding (if on) would be the addition of extra space or time around the data being analysed; in this case, there is none in either direction (unit: seconds)

    cfg.findEvtFree         = 0; % indicates that we are not searching for event-free windows in the data
    cfg.findEvtFreeWin      = 300; % sets the duration of the event-free window to 300

    outSODetect             = slythm.detectSOs_splitStages(cfg, dataEEG); % will be stored in 'outSODetect'

    % Perform spindle detection
    cfg              = [];
    cfg.bpfreq       = [12 18]; % sets a band-pass filter for frequencies between 12 Hz and 18 Hz
    cfg.slStages     = [2 3];

    cfg.thresCh      = dataEEG.label;
    cfg.detectCh     = dataEEG.label;

    cfg.thresType    = 'channelwise';

    cfg.filtType     = 'fir'; % uses a Finite Impulse Response (FIR) filter, another type of filter commonly used in EEG signal processing

    cfg.envelopeType = 'rms'; % applies an RMS (root mean square) envelope to the spindle signal, a common technique to quantify signal amplitude
    cfg.envelopeWin  = 0.2; % window size is set to 0.2 seconds (slides/moves across the signal), smoothing the signal to help detect spindles

    cfg.criterionLen = [0.5, 3]; % spindle duration is constrained to be between 0.5 and 3 seconds

    cfg.criterionCenter     = 'mean';
    cfg.criterionVar        = 'sdCenter'; % normalises the data not just by centring it around its mean, but also scales it according to its standard deviation (divides it by the SD, so the resulting data will have a mean of 0 and a standard deviation of 1)
    cfg.criterionScale      = 1.5;
    cfg.criterionPadding    = [0, 0];

    cfg.searchIndivEvtPeak  = 0; % algorithm will not search for individual peaks within the detected spindle events, instead, it will focus on detecting the spindle events as a whole without dissecting them into smaller peaks
    cfg.searchMargin        = 1; % specifies a time margin (1 second) around detected events to look for additional features or characteristics (algorithm may consider events that occur within 1 second of a detected spindle as relevant)

    cfg.mergeEvts           = 0; % algorithm will not merge events even if they are close in time
    cfg.mergeEvtsGap        = 0.5; % specifies the maximum time gap (in seconds) allowed between events for them to be merged (if this setting is set to 1)

    cfg.doubleThres         = 0; % no double threshold method for event detection (this would involve setting two separate thresholds: one for identifying the beginning of an event and another for its ending)
    cfg.doubleThresScale    = 3; % would indicate that the second threshold is scaled relative to the first threshold, typically set to 3 times the first threshold

    cfg.findEvtFree         = 0;
    cfg.findEvtFreeWin      = 300;

    outfastSpiDetect        = slythm.detectSpindles(cfg, dataEEG);

    % Detect SO-spindle complexes

    SO_Spindle_Cmplx = cell(size(outSODetect.EvtInfo,2), 1); % initialises a cell array to store the detected SO-spindle complexes for each EEG channel. This will have the same number of rows as the number of columns in the outSODetect.EvtInfo structure (which could represent the number of channels or different event types detected), and one column

    % Detect complexes
    for ichann = 1:size(outSODetect.EvtInfo,2) % iterates over each EEG channel
        SO_Spindle_Cmplx_chann = []; % initialises an empty array to hold SO-spindle complex candidates for the current channel
        candidates = outSODetect.EvtInfo(ichann).minTime; % retrieves the minimum time points of detected SO events for the current channel
        for icandi = 1:length(candidates) % begins a nested for loop to iterate over each candidate SO event time for the current channel
            spindle_SO_distances = outfastSpiDetect.EvtInfo(ichann).minTime - candidates(icandi); % calculates the time differences between the detected spindle events and the current SO event candidate (positive values indicate that spindles occur after the SO candidate)

            % Here I check whether the SO has any spindle closer than 1.5
            % seconds from it.
            if any((spindle_SO_distances > 0*new_fs).*(spindle_SO_distances < 1.5*new_fs))
                SO_Spindle_Cmplx_chann = [SO_Spindle_Cmplx_chann, candidates(icandi)]; % appends the current SO candidate to the array (SO_Spindle_Cmplx_chann) for the detected SO-spindle complexes for the current channel
            end
        end
        SO_Spindle_Cmplx{ichann} = SO_Spindle_Cmplx_chann; % saves the collected SO-spindle complexes for the current channel into the cell array
    end

    %% Convert respiration to phase

    % Get respiration frequency of NREM sleep
    mask_NREM = ismember(data_all.scoring{1,1}, [2,3]); % checks which elements in the scoring data are either 2 or 3
    dataResp.trial{1,1} = dataResp.trial{1,1}(:,mask_NREM); % filters the trials to keep only those segments where the subject is in NREM sleep
    dataResp.sampleinfo = [1, sum(mask_NREM)]; % updates the sample information to reflect the new length of the data after filtering (sum of the 1s in the logical array)
    dataResp.time{1,1} = dataResp.time{1,1}(:,mask_NREM); % adjusts the time vector to match the filtered trial data
    
    cfg          = [];
    cfg.length   = 60;
    cfg.overlap  = 0;
    data2_segm    = ft_redefinetrial(cfg, dataResp); % Function from the FieldTrip toolbox to redefine the data trials according to specified parameters

    % run spectral analysis and find peak frequency of respiration
    cfg          = [];
    cfg.channel  = 'Resp';
    cfg.output   = 'pow'; % Requests the power spectrum as the output
    cfg.method   = 'mtmfft'; % Uses multi-taper method for frequency analysis via fast Fourier transform (FFT)
    cfg.taper    = 'hanning'; % Uses a Hanning window for spectral analysis (which helps reduce spectral leakage)
    cfg.foi      = 0.08:0.01:1; % 1/cfg1.length  = 1; % specifies the range of frequencies of interest, from 0.08 Hz to 1 Hz with a step size of 0.01 Hz
    spect        = ft_freqanalysis(cfg, data2_segm); % Performs the frequency analysis on the redefined segments

    [val, ind] = max(spect.powspctrm); % Finds the maximum power value in the spectrum, returning both the value and its index
    peak      = round(spect.freq(ind),2); % Retrieves the frequency corresponding to the maximum power, and rounds this peak frequency to two decimal places
    peak_max  = (peak +0.1); % Define upper frequency bound around the peak frequency (10% higher than the peak frequency)
    peak_min  = (peak -0.1); % Define lower frequency bound around the peak frequency (10% lower than the peak frequency)

    % filter,segement data and derive phase info
    cfg                     = [];
    cfg.detrend             = 'yes'; % Detrends the data to remove linear trends
    cfg.bpfilter            = 'yes'; % Enables bandpass filtering
    cfg.channel             = 'Resp';

    cfg.bpfreq              = [peak_min peak_max]; % Sets the bandpass frequency range based on the peak frequency
    cfg.bpinstabilityfix    = 'reduce'; % Addresses potential instability in the filter design
    cfg.bpfiltdir           = 'twopass';
    data_resp               = ft_preprocessing(cfg, data_all);

    resp_phase              = slythm.doubleInterpolation(data_resp.trial{1,1}(end,:), 10, 0); % Accesses the last trial of the filtered respiration data and interpolates the filtered respiration signal to derive phase information
    respFreqNREM            = peak; % Stores the identified peak frequency of respiration during NREM sleep

    %% Check for circular uniformity test for each channel and store the individual values (CircStat toolbox)
    % For every event, every F (1,2,3,5,z) channel, collect in the table the circular
    % mean and the p-value from the v test
    % for each channel we do a histogramme of the means of each participant
    % once we have that we also do a v test of that (in different
    % sections)

    % Define the list of channels to include
    selectedFrontalChannels = {'F1', 'F2', 'F3', 'F5', 'Fz'};

    % Define event types and their corresponding event times
    % We've done event types at the beginning already
    eventTimes = {outSODetect.EvtInfo, outfastSpiDetect.EvtInfo, SO_Spindle_Cmplx};
    pvals_all = struct();  % Store p-values for all event types

    % Pre-allocate arrays based on column names
    rowValues = cell(1, length(columnNames));

    % Loop through each event type
    for ievent = 1:length(eventTypes)
        event = eventTypes{ievent};
        meanPhases = [];  % Store phase values for averaging
        pValues = [];     % Store p-values for FDR correction
        vValues = [];

        % Select event times based on event type
        switch event
            case 'SO'
                eventTimesForEvent = {outSODetect.EvtInfo.minTime};
            case 'Spindle'
                eventTimesForEvent = {outfastSpiDetect.EvtInfo.minTime};
            case 'SO_Spindle'
                eventTimesForEvent = SO_Spindle_Cmplx;
        end

        % Loop through each channel
        for ichann = 1:length(channels)
            channel = channels{ichann};

            % Check if the current channel is in the list of selected channels
            if ismember(channel, selectedFrontalChannels)

                % Extract respiratory phases for the current event and channel
                eventTimes = eventTimesForEvent{ichann}; % Use the event times from the first channel, since Resp is the only channel
                respPhases = resp_phase(eventTimes);        % Extract respiratory phases using these times

                % Calculate circular mean phase and V-test p-value
                meanPhase = circ_mean(respPhases');
                [pVal, vval] = circ_vtest(respPhases, 0);

                % Store values for later averaging
                meanPhases = [meanPhases, meanPhase];
                pValues = [pValues, pVal];
                vValues = [vValues, vval];

                % Store the results in the pre-allocated arrays
                columnName = sprintf('%s_%s_Mean', eventTypes{ievent}, channels{ichann});
                colIndex = find(strcmp(columnNames, columnName), 1); % Calculate column index
                rowValues{colIndex} = meanPhase;
                rowValues{colIndex + 1} = pVal;
                rowValues{colIndex + 2} = vval;

                % Store p-values in the pvals_all structure for
                % FDR correction
                pvals_all.(event)(ichann) = pVal;
            end
        end

        % Compute the mean of preferred phases across channels
        if ~isempty(meanPhases)
            avgMeanPhase = mean(meanPhases);
        else
            avgMeanPhase = NaN;
        end

        % Apply FDR correction to p-values
        if ~isempty(pValues)
            corrPval_fdr = FDR(pValues', 0.05);  % Indices of significant p-values
            pValues_FDR = pValues;
            pValues_FDR(setdiff(1:length(pValues), corrPval_fdr)) = 0;  % Set non-significant p-values to 0
            avgPVal = mean(pValues);  % Mean of uncorrected p-values
            avgPVal_FDR = mean(pValues_FDR);  % Mean of FDR-corrected p-values
        else
            avgPVal = NaN;
            avgPVal_FDR = NaN;
        end

        rowValues{find(strcmp(columnNames, [event '_Mean_Phase']), 1)} = avgMeanPhase;
        rowValues{find(strcmp(columnNames, [event '_Mean_pVal']), 1)} = avgPVal;
        rowValues{find(strcmp(columnNames, [event '_Mean_pVal_FDR']), 1)} = avgPVal_FDR;

    end

    % Store the computed values in the table
    rowValues{find(strcmp(columnNames, 'SubjectID'), 1)} = subjectID;
    rowValues{find(strcmp(columnNames, 'Age'), 1)} = 25;
    rowValues{find(strcmp(columnNames, 'Gender'), 1)} = 'F';

    %%

    % Create a table for storing results of all subjects
    newRow = array2table(rowValues, 'VariableNames', columnNames);

    newRow.Wake = stageDurations.Wake;
    newRow.N1 = stageDurations.N1;
    newRow.N2 = stageDurations.N2;
    newRow.N3 = stageDurations.N3;
    newRow.REM = stageDurations.REM;

    % Calculate total sleep time
    newRow.TST = stageDurations.N1 + stageDurations.N2 + stageDurations.N3 + stageDurations.REM;

    subjectsTable = [subjectsTable; newRow];  % Append newRow to the pre-defined subjectsTable

end

%% Save the overall table to a .mat file

save('subjectsTable.mat', 'subjectsTable');
writetable(subjectsTable, 'subjectsTable.xlsx');