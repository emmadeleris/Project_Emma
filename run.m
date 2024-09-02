%% Set paths and code variables
clear all

subjectsFolder = "";
subjects = dir(subjectsFolder);

%% For loop

for isubject=1:height(subjects)

    %% Load the data (dont forget the hypno!)
    data = edfread(); % Try to read header too somehow

    %% See if there is a way to convert basic edf to ft structure

    %% Preprocessing (e.g. changing reference to M1 and M2)
    % VERY IMPORTANT: Preprocess EEG and Respiration separately

    %% Detect events IN THE EEG (specific channels) - check detecting algos from Esteban

    %% Convert respiration to phase

    %% Check the phases where the SO and Spindles are occuring

    %% Check for circular means and circular uniformity test (CircStat toolbox)

    %% We store relevant information from EACH subject

end

%% Check population results
% Circ plot the means from all subjects

%% Run population stats
% Another Rayleight test 