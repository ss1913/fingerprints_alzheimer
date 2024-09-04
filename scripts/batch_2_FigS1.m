%% Batch for Figure 2 and Supplementary Figure 1
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

% Load data for GENEVA
FC_GENEVA = load(fullfile(Path2Inputs, 'FCs_GENEVA.mat'));
T_GENEVA = readtable(fullfile(Path2Inputs, 'Behav_GENEVA.csv'));
Motion_GENEVA = load(fullfile(Path2Inputs, 'Motion_GENEVA.mat'));

% Load data for ADNI
FC_ADNI = load(fullfile(Path2Inputs, 'FCs_ADNI.mat'));
T_ADNI = readtable(fullfile(Path2Inputs, 'Behav_ADNI.csv'));
Motion_ADNI = load(fullfile(Path2Inputs, 'Motion_ADNI.mat'));

% Load parcellation
 parc = load(fullfile(Path2Inputs, 'shen_yeo_RS7.mat'));

% Load results
load(fullfile(Path2Outputs, 'ID_metrics.mat'));

%% Define constants and variables

% Constants
N = size(parc.yeoROIs, 1);
groups = {'CU_minus', 'MCI_plus', 'Dementia_plus'};
groups_labels = {'CU AB-', 'MCI AB+', 'AD Dementia'};
datasets = {'GENEVA', 'ADNI'};
% Masks
mask_ut = triu(true(N), 1); % Upper triangle mask
cerebROI = parc.yeoROIs == 9; % ROIs for Cerebellum 
mask_cereb = bsxfun(@or, cerebROI, cerebROI.'); % Mask for Cerebellum
mask_cereb_ut = mask_cereb(mask_ut); % Upper triangle mask for Cerebellum


%% Save data for statistics in R

for j = 1:length(datasets)
    dataset_name = datasets{j};
    Behav = eval(['T_', dataset_name]); 
    Motion = eval(['Motion_', dataset_name]);   
    % Initialize arrays for concatenation
    g = [];
    iself = [];
    iothers = [];
    id = [];
    age = [];
    sex = [];
    yoe = [];
    scanner = [];
    delta_fd = [];
    fd = [];    
    for i = 1:length(groups)
        group_name = groups{i};
        % Create vector for groups
        g = [g, repmat({group_name}, 1, length(ID_metrics.(dataset_name).(group_name).Iself_vect))];
        % Create vector for Iself
        iself = [iself; ID_metrics.(dataset_name).(group_name).Iself_vect];
        % Create vector for IOthers
        iothers = [iothers; ID_metrics.(dataset_name).(group_name).Iothers_vect];        
        % Create a vector for ID
        id = [id; Behav.ID(strcmp(Behav.Group, group_name))];
        % Create a vector for Age
        age = [age; Behav.Age(strcmp(Behav.Group, group_name))];       
        % Create a vector for Sex
        sex = [sex; Behav.Sex(strcmp(Behav.Group, group_name))];
        % Create a vector for YoE
        yoe = [yoe; Behav.PTEDUCAT(strcmp(Behav.Group, group_name))];      
        % Create a vector for Scanner type if the dataset is ADNI
        if strcmp(dataset_name, 'ADNI')
            scanner = [scanner; Behav.Machine(strcmp(Behav.Group, group_name))];
        end   
        % Create a vector for head motion (FD)     
        fd = [fd; Motion.FD(strcmp(Motion.Group,group_name))];                
        % Create a vector for delta_FD
        temp = Motion.FD_test - Motion.FD_retest;
        delta_fd = [delta_fd; temp(strcmp(Motion.Group,group_name))];
    end
    % Store the concatenated vector in a structure
    group_struct.(dataset_name) = g';
    iself_struct.(dataset_name) = iself;
    iothers_struct.(dataset_name) = iothers;
    id_struct.(dataset_name) = id;
    age_struct.(dataset_name) = age;
    sex_struct.(dataset_name) = sex;
    yoe_struct.(dataset_name) = yoe;
    fd_struct.(dataset_name) = fd;
    delta_fd_struct.(dataset_name) = delta_fd;
    if strcmp(dataset_name, 'ADNI')
        scanner_struct.(dataset_name) = scanner;
    end
    % Create table for Iothers
    if strcmp(dataset_name, 'ADNI')
        tbl_iothers.(dataset_name) = table(iothers_struct.(dataset_name), age_struct.(dataset_name), sex_struct.(dataset_name), yoe_struct.(dataset_name), fd_struct.(dataset_name), scanner_struct.(dataset_name), group_struct.(dataset_name), id_struct.(dataset_name), ...
            'VariableNames', {'Iothers', 'Age', 'Sex', 'YoE', 'FD', 'Scanner', 'Group', 'Subject'});
    else
        tbl_iothers.(dataset_name) = table(iothers_struct.(dataset_name), age_struct.(dataset_name), sex_struct.(dataset_name), yoe_struct.(dataset_name), fd_struct.(dataset_name),group_struct.(dataset_name), id_struct.(dataset_name), ...
            'VariableNames', {'Iothers', 'Age', 'Sex', 'YoE', 'FD','Group', 'Subject'});
    end
    writetable(tbl_iothers.(dataset_name), fullfile(Path2Outputs, ['Iothers_' dataset_name '.csv']));
    % Create table for Iself
    if strcmp(dataset_name, 'ADNI')
        tbl_iself.(dataset_name) = table(iself_struct.(dataset_name), age_struct.(dataset_name), sex_struct.(dataset_name), yoe_struct.(dataset_name), delta_fd_struct.(dataset_name), scanner_struct.(dataset_name), group_struct.(dataset_name), id_struct.(dataset_name), ...
            'VariableNames', {'Iself', 'Age', 'Sex', 'YoE', 'delta_FD','Scanner', 'Group', 'Subject'});
    else
        tbl_iself.(dataset_name) = table(iself_struct.(dataset_name), age_struct.(dataset_name), sex_struct.(dataset_name), yoe_struct.(dataset_name), delta_fd_struct.(dataset_name), group_struct.(dataset_name), id_struct.(dataset_name), ...
            'VariableNames', {'Iself', 'Age', 'Sex', 'YoE', 'delta_FD','Group', 'Subject'});
    end
    writetable(tbl_iself.(dataset_name), fullfile(Path2Outputs, ['Iself_' dataset_name '.csv']));
end

%% Supplementary Figure 1

colors = lines(length(groups)); % Generate different colors for each group

for j = 1:length(datasets)
    dataset_name = datasets{j};
    figure;
    % Define data sources and labels
    tables = {'tbl_iself', 'tbl_iothers'};
    data_sources = {'Iself', 'Iothers'};
    y_labels = {'iSelf', 'iOthers'};

    for k = 1:2
        subplot (1,2,k);
        data_source = data_sources{k};
        y_label = y_labels{k};
        table = eval(tables{k});
        % Extract data for current dataset
        data = table.(dataset_name).(data_source);
        group = table.(dataset_name).Group;
        % Create boxplot
        boxplot(data, group, 'Notch', 'on');
        hold on;
        % Add individual data points
        for i = 1:length(groups)
            current_group_data = data(strcmp(group, groups{i}));
            scatter(repmat(i, length(current_group_data), 1), current_group_data, ...
                'filled', 'MarkerFaceColor', colors(i,:));
        end
        % Add labels and title
        xlabel('Group');
        ylabel(y_label);
        title(dataset_name);
        axis square
        % Add specific caption for each subplot
        if k == 1
            title(['A) ', dataset_name]);
        else
            title(['B) ', dataset_name]);
        end
        hold off;
        sgtitle('Supplementary Figure 1');
    end
end

% see batch_3_TableS1.R for one-way ANOVAs to test the effect of group on ISelf and IOthers separately,
% after checking for nuisance variables  (Age, Sex, YoE and absolute difference between motion (FD) at test vs. retest), 
% with 5000 permutations to control for sample size differences.
