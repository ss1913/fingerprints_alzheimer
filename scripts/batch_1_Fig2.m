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

%% Define constants and variables

% Constants
N = size(parc.yeoROIs, 1);
null_runs = 1000; 
groups = {'CU_minus', 'MCI_plus', 'Dementia_plus'};
groups_labels = {'CU AB-', 'MCI AB+', 'AD Dementia'};
datasets = {'GENEVA', 'ADNI'};
% Masks
mask_ut = triu(true(N), 1); % Upper triangle mask
cerebROI = parc.yeoROIs == 9; % ROIs for Cerebellum 
mask_cereb = bsxfun(@or, cerebROI, cerebROI.'); % Mask for Cerebellum
mask_cereb_ut = mask_cereb(mask_ut); % Upper triangle mask for Cerebellum

%% Compute head motion statistics

motion_stats = struct();
for j = 1:length(datasets)
    dataset_name = datasets{j};
    Motion = eval(['Motion_', dataset_name]); 

    for i = 1:length(groups)
        group_name = groups{i};
        g = zeros(size(Motion.Group, 1), 1); 

         % Create a numerical vector for the groups
         for i = 1:length(groups)
             group_name = groups{i};
             g(strcmp(Motion.Group, group_name)) = i;
         end
    end
    % Head motion across groups
    motion_stats.(dataset_name).FD_acrossgroups_pval = anova1(Motion.FD, g);
    motion_stats.(dataset_name).TaggedVols_acrossgroups_pval = anova1(Motion.perc_tag_vols, g);
    % Head motion across test and retest
    [~, motion_stats.(dataset_name).FD_test_retest_pval] = ttest(Motion.FD_test, Motion.FD_retest);
    [~, motion_stats.(dataset_name).TaggedVols_test_retest_pval] = ttest(Motion.perc_tag_vols_test, Motion.perc_tagvols_retest);

end

close all

%% Compute ID matrices and metrics (ISelf, IOthers, IDiff, Success-rate)

%ID matrices
ID_mats = struct();
for j = 1:length(datasets)
    dataset_name = datasets{j};
    FC = eval(['FC_', dataset_name]);    
    mask_nan_ut = FC.mask_nan_90(mask_ut); % Upper triangle NaN mask
    valid_mask = ~(mask_nan_ut | mask_cereb_ut); % Valid data points mask
    ID_mats.(dataset_name) = struct();    
    for i = 1:length(groups)
        group_name = groups{i};
        test = ['FC_2D_test_', groups{i}];
        retest = ['FC_2D_retest_', groups{i}];
        ID_mats.(dataset_name).(group_name) = corr(FC.(test)(valid_mask,:), FC.(retest)(valid_mask,:), 'rows', 'pairwise');
    end
end

% ID metrics
ID_metrics_names = {'idiff_m', 'idiff_se', 'idiff_vect', 'success', 'iself_m', 'iself_sd','Iself_vect', 'Iothers_vect', 'Iothers_m'};
ID_metrics = struct();
ID_metrics_null = struct();

for j = 1:length(datasets)
    dataset_name = datasets{j};
    ID_metrics.(dataset_name) = struct(); % Initialize as a structure
    for i = 1:length(groups)
        group_name = groups{i};
        [ID_metrics.(dataset_name).(group_name).idiff_m, ...
         ID_metrics.(dataset_name).(group_name).idiff_se, ...
         ID_metrics.(dataset_name).(group_name).idiff_vect, ...
         ID_metrics.(dataset_name).(group_name).success, ...
         ID_metrics.(dataset_name).(group_name).iself_m, ...
         ID_metrics.(dataset_name).(group_name).iself_sd, ...
         ID_metrics.(dataset_name).(group_name).Iself_vect, ...
         ID_metrics.(dataset_name).(group_name).Iothers_vect, ...
         ID_metrics.(dataset_name).(group_name).Iothers_m] = f_ID(ID_mats.(dataset_name).(group_name));
    end
    ID_metrics_null.(dataset_name) = struct(); % Initialize as a structure
    for i = 1:length(groups)
        group_name = groups{i};
        [ID_metrics_null.(dataset_name).(group_name).idiff_null, ...
         ID_metrics_null.(dataset_name).(group_name).success_null] = f_ID_permutation(ID_mats.(dataset_name).(group_name), null_runs);
        % Test if real IDiff and Success-rate significantly different than chance at 0.01 
        thr = prctile(ID_metrics_null.(dataset_name).(group_name).idiff_null, 99);
        results.(dataset_name).(group_name).idiff_sig = (ID_metrics.(dataset_name).(group_name).idiff_m>thr);
        thr = prctile(ID_metrics_null.(dataset_name).(group_name).success_null, 99);
        results.(dataset_name).(group_name).success_sig = (ID_metrics.(dataset_name).(group_name).success>thr);
    end
end

save(fullfile(Path2Outputs, 'ID_metrics.mat'));

%% Figure 2 

resolution = 600;

for j = 1:length(datasets)
    dataset_name = datasets{j};
    
    % ID matrices - Figures 2A and 2C
    figure;
    for i = 1:length(groups)
        group_name = groups{i}; 
        subplot(1, 3, i);
        imagesc(ID_mats.(dataset_name).(group_name), [0.2, 0.7]);
        axis square;
        title(['ID ', groups_labels{i}]);
        xlabel('Subj test');
        ylabel('Subj retest');
        colorbar;
    end
    if strcmp(dataset_name, 'GENEVA')
        sgtitle(['Figure 2A - ', 'GENEVA']);
    else
        sgtitle(['Figure 2C - ', 'ADNI']);
    end
   

    % Boxplots - Figures 2B and 2D
    figure;
    for i = 1:length(groups)
        group_name = groups{i};
        subplot(1, 3, i);
        measures = cat(2, ID_metrics.(dataset_name).(group_name).Iself_vect, ID_metrics.(dataset_name).(group_name).Iothers_vect);
        boxplot(measures, 'Symbol', 'k.'); hold on;
        title([groups_labels{i}]);
        xticklabels({'ISelf', 'IOthers'});
        ylim([0.15, 1]);
        parallelcoords(measures, 'Color', 0.8 * [1, 1, 1], 'LineStyle', '-','Marker', '.', 'MarkerSize', 10);
        axis square;
    end
    if strcmp(dataset_name, 'GENEVA')
        sgtitle(['Figure 2B - ', 'GENEVA']);
    else
        sgtitle(['Figure 2D - ', 'ADNI']);
    end
    

end





