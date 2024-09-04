%% Batch for Figure X
%Author: Sara Stampacchia
% 
clear
clc
close all

%% Initialize paths and load data

% Define paths
Path2Inputs = fullfile('..','inputs/');
Path2Outputs = fullfile('..','outputs/');
addpath(fullfile('..', 'functions/'));

% Load parcellation
parc = load(fullfile(Path2Inputs, 'shen_yeo_RS7.mat'));

% Load results
load(fullfile(Path2Outputs, "ICC_mats.mat"));


%% Define constants and variables
% Constants
groups = {'CU_minus', 'MCI_plus', 'Dementia_plus'};
datasets = {'GENEVA', 'ADNI'};


%% Supplementary Figure 4

%Compute nodal strenght for brain renders 
for j = 1:length(datasets)
    dataset_name = datasets{j};
    for i = 1:length(groups)
            group_name = groups{i};
            ICC_mats.(dataset_name).(group_name).nodalstrenght = f_plot_Bnet_renders(ICC_mats.(dataset_name).(group_name).ICC_mat_sig, Path2Outputs, ['ICC_sig_', dataset_name, '_' ,group_name, '_'],1,Path2Inputs);
            ICC_mats.(dataset_name).(group_name).nodalstrenght_percentiles_75_100 = prctile(ICC_mats.(dataset_name).(group_name).nodalstrenght,[75 100]);
    end
end

save(fullfile(Path2Outputs, 'ICC_mats.mat'), 'ICC_mats', '-v7.3');

% Brain render were visualised using BrainNet viewer - options for loading
% are provided in Path2Inputs, 'options_brain_renders.mat'. 
% For visualisation input the value corresponding to the 75th percentile for each group
% in the relevant option field.



