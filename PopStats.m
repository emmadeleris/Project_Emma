%% Set paths and code variables
% This needs to go directly to the subjectsTable here – to modify
clear all

populationFolder = "/Users/emma/Documents/GitHub/Project_Emma";
load subjectsTable.mat

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
        ColumnNames{end+1} = [eventTypes{ievent}, '_', channels{ichann}, '_Plot'];
    end
end

% Create the population table
populationTable = cell2table(cell(0, length(ColumnNames)), 'VariableNames', ColumnNames);

% Loop through event types and channels to calculate across-participant stats
rowValues = cell(1, length(ColumnNames)) % rowValues = NaN(1, length(ColumnNames)); % Initialise a row for results
for ievent = 1:length(eventTypes)
    event = eventTypes{ievent};
    for ichann = 1:length(channels)
        channel = channels{ichann};

        % Construct column names for mean phases across participants
        meanColumn = [event, '_', channel, '_Mean'];

        % Extract mean phases across participants
        participantMeans = subjectsTable.(meanColumn); % Get participant means
        participantMeans = participantMeans(~isnan(participantMeans)); % Exclude NaNs

        % Calculate the circular mean
        popCircMean = circ_mean(participantMeans);

        % Perform Rayleigh test
        [popPVal, ~] = circ_rtest(participantMeans);

         % Create the circular plot
        fig = figure('Visible', 'off'); % Generate figure but keep it hidden (can also be changed)
        circ_plot(participantMeans, 'hist', [], 20, true, true, 'linewidth', 2);
        title([event, ' - ', channel]);

        % Store the results
        colIndexMean = find(strcmp(ColumnNames, meanColumn));
        colIndexPVal = colIndexMean + 1; % pVal column is next to the Mean column
        colIndexPlot = colIndexMean + 2; % Plot column is next to pVal column
        rowValues{colIndexMean} = popCircMean;
        rowValues{colIndexPVal} = popPVal;
        rowValues{colIndexPlot} = fig;
    end
end

% Add the calculated row to the summary table
popRow = cell2table(rowValues, 'VariableNames', ColumnNames); % summaryRow = array2table(rowValues, 'VariableNames', ColumnNames);
populationTable = [populationTable; popRow];

% Display or save the summary table
save('populationTable.mat','populationTable');
writetable(populationTable, 'populationTable.xlsx')

% Display a figure
figure(populationTable{1, 'SO_F1_Plot'}); % Amend line as needed

%% Run population stats
% Run statistical analysis on the population-level results (e.g., Rayleigh test across subjects)
% Perform Rayleigh test for circular uniformity on the population means

%% Check population results
% Circ plot the means from all subjects

%% Run population stats
% Another Rayleigh test