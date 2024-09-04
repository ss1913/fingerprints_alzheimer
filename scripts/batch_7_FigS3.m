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

% Load data for GENEVA
FC_GENEVA = load(fullfile(Path2Inputs, 'FCs_GENEVA.mat'));
T_GENEVA = readtable(fullfile(Path2Inputs, 'Behav_GENEVA.csv'));

% Load data for ADNI
FC_ADNI = load(fullfile(Path2Inputs, 'FCs_ADNI.mat'));
T_ADNI = readtable(fullfile(Path2Inputs, 'Behav_ADNI.csv'));

% Load parcellation
parc = load(fullfile(Path2Inputs, 'shen_yeo_RS7.mat'));

% Load results
load(fullfile(Path2Outputs, "ICC_mats.mat"));

%% Define constants and variables
% Constants
N = size(parc.yeoROIs, 1);
num_runs = 10; % To reduce computational time, here num_runs is here set to 10. To replicate manuscript set to 1000 (~8hrs of computational time).
groups = {'CU_minus', 'MCI_plus', 'Dementia_plus'};
groups_labels = {'CU AB-', 'MCI AB+', 'AD Dementia'};
datasets = {'GENEVA', 'ADNI'};
numedges = size(FC_GENEVA.FC_2D_retest_CU_minus,1);
% Masks
mask_ut = triu(true(N), 1); % Upper triangle mask
cerebROI = parc.yeoROIs == 9; % ROIs for Cerebellum 
mask_cereb = bsxfun(@or, cerebROI, cerebROI.'); % Mask for Cerebellum
mask_cereb_ut = mask_cereb(mask_ut); % Upper triangle mask for Cerebellum



%% Compute ICC null

% --> Alert: To reduce computational time, here num_runs is here set to 10 (see line 30). 
% To replicate manuscript results set num_runs to 1000 (this will take ~8hrs). 
% Results with 1000 runs are saved in outputs/manuscript

ICC_mats_null = struct();
for j = 1:length(datasets)
    tic
    dataset_name = datasets{j};
    disp(dataset_name)
    FC = eval(['FC_', dataset_name]);
    mask_nan_ut = FC.mask_nan_90(mask_ut); % Upper triangle NaN mask
    valid_mask = (FC.mask_nan_90|mask_cereb); % remove cerebellum and mask nans
    if strcmp(dataset_name, 'GENEVA')
        n_subset = 5;
    else
        n_subset = 10;
    end
    ICC_mat_3D_null = nan(N,N,num_runs);
    for run=1:num_runs
        disp(['run ', num2str(run)]);
        aux_test_null_groups = zeros(numedges,n_subset,length(groups));
        aux_retest_null_groups = zeros(numedges,n_subset,length(groups));    
        for i = 1:length(groups)
            group_name = groups{i};
            test = ['FC_2D_test_', groups{i}];
            retest = ['FC_2D_retest_', groups{i}];
            % Lottery
            N_sub = size(FC.(test),2);
            lottery = randperm(N_sub,n_subset);
            % Select random subjects from FC test and retest from group
            aux_test_null = FC.(test)(:,lottery);
            aux_retest_null= FC.(retest)(:,lottery);
            % Iteratively save FC test and retest for current group
            aux_test_null_groups(:,:,i) = aux_test_null;
            aux_retest_null_groups(:,:,i) = aux_retest_null;    
        end
        % Create a null FC including random subjects from each group
        test_null = reshape(aux_test_null_groups, numedges, []);
        retest_null = reshape(aux_retest_null_groups, numedges, []);
        % Compute ICC null
        [~, ICC_mat_null] = f_ICC_edgewise(test_null,retest_null,N);
        ICC_mat_null(valid_mask) = 0; % remove cerebellum and mask nans
        ICC_mat_3D_null(:,:,run) = ICC_mat_null;
    end
    ICC_mats_null.(dataset_name).ICC_mat_null_3D = ICC_mat_3D_null;
    ICC_mats_null.(dataset_name).ICC_mat_null = nanmean(ICC_mats_null.(dataset_name).ICC_mat_null_3D,3);
end
disp('Done!');
toc

save(fullfile(Path2Outputs, 'ICC_mats_null.mat'), 'ICC_mats_null', '-v7.3');


%% Compute edge-wise p-value
%Significant edges: p-value represents proportion of times across runs when 
%edgewise-null is bigger than edgewise-real (i.e., probabily of accepting the null hypothesis)

for j = 1:length(datasets)
    dataset_name = datasets{j};
    for i = 1:length(groups)
        group_name = groups{i};
        pos = zeros(N,N, num_runs);
        for run=1:num_runs
            aux = ICC_mats_null.(dataset_name).ICC_mat_null_3D(:,:,run);
            pos(:,:,run) = ICC_mats.(dataset_name).(group_name).ICC_mat > aux;
        end
        p_mat = 1-(sum(pos,3)/num_runs);
        p_mat_sig = p_mat < 0.05;
        sum(p_mat_sig(mask_ut))/numedges*100
        ICC_mats.(dataset_name).(group_name).ICC_mat_sig = ICC_mats.(dataset_name).(group_name).ICC_mat.*p_mat_sig;
    end
end

save(fullfile(Path2Outputs, 'ICC_mats.mat'), 'ICC_mats', '-v7.3');

%% Supplementary Figure 3 
% Define ordered Yeo ROIs and exclude cerebellum ROIs (cerebellum = 9)
yeoROIs_ordered = [(1:N)', (parc.yeoROIs(parc.yeoOrder))]; 
index = (find(yeoROIs_ordered(:, 2) == 9, 1))-1;
ROIs = (1:index); % To remove  cerebellum 

range = [0.4, 0.8];

for j = 1:length(datasets)
    dataset_name = datasets{j};
    % ID matrices - Figures 2A and 2C
    figure;
    hold on
    subplot(1, 4, 1);
    imagesc(ICC_mats_null.(dataset_name).ICC_mat_null(parc.yeoOrder(ROIs), parc.yeoOrder(ROIs)), range);
        axis square;
        xlabel('ROI test');
        ylabel('ROI retest');
        title(['A) ICC Surrogate group unspecific', groups_labels{i}]);
        colorbar;
                
    for i = 1:length(groups)
        group_name = groups{i}; 
        subplot(1, 4, i+1);
        imagesc(ICC_mats.(dataset_name).(group_name).ICC_mat_sig(parc.yeoOrder(ROIs),parc.yeoOrder(ROIs)) ,range);
        axis square;
        xlabel('ROI test');
        ylabel('ROI retest');
        title(['B) Significant ICC matrices ', groups_labels{i}]);
        colorbar;
    end
    if strcmp(dataset_name, 'GENEVA')
        sgtitle(['Supplementary Figure 3 - ', 'GENEVA']);
    else
        sgtitle(['Supplementary Figure 3 - ', 'ADNI']);
    end
end
 


