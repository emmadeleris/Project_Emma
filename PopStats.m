%% Set paths and code variables
clear all

populationFolder = "/Users/emma/Documents/GitHub/Project_Emma";
load subjectsTable.mat

base_dir = "/Users/emma/Desktop/München/MSc Learning Sciences/RP Schreiner/Data analysis";
fieldtrip_path = "toolbox/fieldtrip-20240916";
breathmetrics_path = "toolbox/breathmetrics";
circstat_path = "toolbox/CircStat2012a";

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

%% Check population results
% After processing all subjects, we now check and analyse the population results
% for circular mean and p-values across subjects.

% Initialise a population-level table
populationTable = table();

% Extract unique event types and channels
eventTypes = {'SO', 'Spindle', 'SO_Spindle'};
channels = {'F1', 'F2', 'F3', 'F5', 'Fz'};

% Predefine column names for the population table
ColumnNames = {};
for ievent = 1:length(eventTypes)
    for ichann = 1:length(channels)
        ColumnNames{end+1} = [eventTypes{ievent}, '_', channels{ichann}, '_Mean'];
        ColumnNames{end+1} = [eventTypes{ievent}, '_', channels{ichann}, '_pVal'];
        ColumnNames{end+1} = [eventTypes{ievent}, '_', channels{ichann}, '_VVal']; % Vector strength
        ColumnNames{end+1} = [eventTypes{ievent}, '_', channels{ichann}, '_CircStd']; % Circular standard deviation
        ColumnNames{end+1} = [eventTypes{ievent}, '_', channels{ichann}, '_Plot'];
    end
end

% Add mean phase, mean p-value, and vector strength columns for each event type
for ievent = 1:length(eventTypes)
    ColumnNames{end+1} = [eventTypes{ievent}, '_Mean_Phase'];
    ColumnNames{end+1} = [eventTypes{ievent}, '_Mean_pVal'];
    ColumnNames{end+1} = [eventTypes{ievent}, '_Mean_VVal']; % Mean vector strength for cross-electrode
    ColumnNames{end+1} = [eventTypes{ievent}, '_Mean_pVal_FDR']; % FDR-corrected p-value
    ColumnNames{end+1} = [eventTypes{ievent}, '_Mean_CircStd']; % Standard deviation for cross-electrode
    ColumnNames{end+1} = [eventTypes{ievent}, '_All_Electrodes_Plot']; % Plot for cross-electrode
end

% Create the population table
populationTable = cell2table(cell(0, length(ColumnNames)), 'VariableNames', ColumnNames);

% Loop through event types and channels to calculate across-participant stats
rowValues = cell(1, length(ColumnNames)); % Initialise a row for results
for ievent = 1:length(eventTypes)
    event = eventTypes{ievent};
    for ichann = 1:length(channels)
        channel = channels{ichann};

        % Construct column names for mean phases across participants
        meanColumn = [event, '_', channel, '_Mean'];
        overallMeanColumn = [event, '_Mean_Phase'];

        % Extract mean phases across participants
        participantMeans = [subjectsTable.(meanColumn){:}]; % Get participant means
        participantOverallMeans = [subjectsTable.(overallMeanColumn){:}]; % Get participant cross-electrode means
        participantMeans = participantMeans(~isnan(participantMeans)); % Exclude NaNs
        participantOverallMeans = participantOverallMeans(~isnan(participantOverallMeans)); % Exclude NaNs

        % Calculate the circular mean and circular standard deviation
        popCircMean = circ_mean(participantMeans');
        popCircStdv = circ_std(participantMeans');
        popOverallCircMean = circ_mean(participantOverallMeans');
        popOverallCircStdv = circ_std(participantOverallMeans');

        % Perform v test
        [popPVal, popVVal] = circ_vtest(participantMeans, 0);
        [OverallPopPVal, OverallPopVVal] = circ_vtest(participantOverallMeans, 0);

        % Create the circular plot for single-electrode data
        SingleElectFigure = figure();
        h = circ_plot(participantMeans', 'hist', [], 20, true, true, 'linewidth', 4, 'color', 'k');
        hlines = findall(gcf, 'Type', 'line');
        set(hlines(2), 'LineWidth', 5, 'Color', 'r');
        set(hlines(3:11), 'LineWidth', 1);
        title(strrep([event, ' - ', channel], '_', '-'));

        % Create the circular plot for cross-electrode data
        CrossElectFigure = figure();
        h = circ_plot(participantOverallMeans', 'hist', [], 20, true, true, 'linewidth', 4, 'color', 'k');
        hlines = findall(gcf, 'Type', 'line');
        set(hlines(2), 'LineWidth', 5, 'Color', 'r');
        set(hlines(3:11), 'LineWidth', 1);
        title(strrep([event, ' - All Electrodes'], '_', '-'));

        % Save the circular plot in a dedicated folder
        figDir = "/Users/emma/Documents/GitHub/Project_Emma/Population figures";

        % Saving the single-electrode figure
        SingleElectFigname = fullfile(figDir, [event, '_', channel]);
        saveas(SingleElectFigure, SingleElectFigname, 'jpeg');

        % Saving the cross-electrode figure
        CrossElectFigname = fullfile(figDir, [event, '_All_Electrodes']);
        saveas(CrossElectFigure, CrossElectFigname, 'jpeg');

        % Store the results for the single-electrode
        colIndexSingleElectMean = find(strcmp(ColumnNames, meanColumn));
        colIndexSingleElectPVal = colIndexSingleElectMean + 1;
        colIndexSingleElectVVal = colIndexSingleElectMean + 2;
        colIndexSingleElectStdv = colIndexSingleElectMean + 3;
        colIndexSingleElectPlot = colIndexSingleElectMean + 4;

        rowValues{colIndexSingleElectMean} = popCircMean;
        rowValues{colIndexSingleElectPVal} = popPVal;
        rowValues{colIndexSingleElectVVal} = popVVal; % Store vector strength
        rowValues{colIndexSingleElectStdv} = popCircStdv; % Store circular standard deviation
        rowValues{colIndexSingleElectPlot} = SingleElectFigure;

        % Store the results for the cross-electrode
        colIndexCrossElectMean = find(strcmp(ColumnNames, [event, '_Mean_Phase']));
        colIndexCrossElectPVal = colIndexCrossElectMean + 1;
        colIndexCrossElectVVal = colIndexCrossElectMean + 2;
        colIndexCrossElectFDR = colIndexCrossElectMean + 3;
        colIndexCrossElectStdv = colIndexCrossElectMean + 4;
        colIndexCrossElectPlot = colIndexCrossElectMean + 5;

        rowValues{colIndexCrossElectMean} = popOverallCircMean;
        rowValues{colIndexCrossElectPVal} = OverallPopPVal;
        rowValues{colIndexCrossElectVVal} = OverallPopVVal; % Store vector strength for cross-electrode
        rowValues{colIndexCrossElectStdv} = popOverallCircStdv; % Store circular standard deviation for cross-electrode
        rowValues{colIndexCrossElectPlot} = CrossElectFigure;
    end % End of inner loop (channels)

    %% Apply FDR correction to p-values
    % Assign pValues to either popPVal or OverallPopPVal (this was actually
    % silly because I'm only working with OverallPopPVal oops)
    pValues = popPVal;  % Default to popPVal
    if exist('OverallPopPVal', 'var') && ~isempty(OverallPopPVal)
        pValues = OverallPopPVal;  % Use OverallPopPVal if available
    end

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

    % Store the FDR-corrected p-values for the cross-electrode analyses
    colIndexCrossElectFDR = find(strcmp(ColumnNames, [event, '_Mean_pVal_FDR']));
    rowValues{colIndexCrossElectFDR} = avgPVal_FDR;
end % End of outer loop (eventTypes)

% Add the calculated row to the summary table
popRow = cell2table(rowValues, 'VariableNames', ColumnNames);
populationTable = [populationTable; popRow];

% Save the population table
save('populationTable.mat', 'populationTable');
writetable(populationTable, 'populationTable.xlsx');
