%% Batch for Figure 3
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

%% Compute ICC

% --> Alert: To reduce computational time, here num_runs is here set to 10 (see line 30). 
% To replicate manuscript results set num_runs to 1000 (this will take ~8hrs). 
% Results with 1000 runs are saved in outputs/manuscript

ICC_mats = struct();
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
    
    for i = 1:length(groups)
        group_name = groups{i};
        disp(group_name);
        ICC_mat_3D = nan(N,N,num_runs);
        for run=1:num_runs
            disp(['run ', num2str(run)]);
            test = ['FC_2D_test_', groups{i}];
            retest = ['FC_2D_retest_', groups{i}];
            lottery = randperm(size(FC.(test),2),n_subset);
            [~, ICC_mat] = f_ICC_edgewise(FC.(test)(:,lottery),FC.(retest)(:,lottery),N);
            ICC_mat(valid_mask) = 0; % remove cerebellum and mask nans
            ICC_mats.(dataset_name).(group_name).ICC_mat_3D(:,:,run) = ICC_mat; % save it in the structure
        end
        
        % average across bootstrap runs
        aux = mean(ICC_mats.(dataset_name).(group_name).ICC_mat_3D,3);
        aux(aux==0)=nan; % set 0 to NaN so that it does influence Min value as well as percentile calculation
        ICC_mats.(dataset_name).(group_name).ICC_mat = aux;
    end
    toc
end
disp('Done!')


%% Create overlap GENEVA + ADNI

thr = 0.6; %threshold at ICC=0.6
for j = 1:length(datasets)
    dataset_name = datasets{j};
    for i = 1:length(groups)
        group_name = groups{i};
        % Binarize ICC matrices
        ICC_mat_thr = zeros(N,N);
        ICC_mat_thr(ICC_mats.(dataset_name).(group_name).ICC_mat>thr)=1;
        ICC_mats.(dataset_name).(group_name).ICC_mat_thr=ICC_mat_thr; % save it in the structure 
    end
end

for i = 1:length(groups)
    group_name = groups{i};
    % Create overlap 
    overlap = ICC_mats.GENEVA.(group_name).ICC_mat_thr + ICC_mats.ADNI.(group_name).ICC_mat_thr;
    overlap = overlap>=2;
    ICC_mats.GENEVA_ADNI.(group_name).overlap = overlap;
end

save(fullfile(Path2Outputs, 'ICC_mats.mat'), 'ICC_mats', '-v7.3');



%% Figure 3

% Define ordered Yeo ROIs and exclude cerebellum ROIs (cerebellum = 9)
yeoROIs_ordered = [(1:N)', (parc.yeoROIs(parc.yeoOrder))]; 
index = (find(yeoROIs_ordered(:, 2) == 9, 1))-1;
ROIs = (1:index); % To remove  cerebellum 

% Loop through each dataset and group to create figures 3A and 3B
for j = 1:length(datasets)
    dataset_name = datasets{j};
    figure;

    for i = 1:length(groups)
        group_name = groups{i};
        % Calculate percentile range for ICC matrix
        range = [prctile(ICC_mats.(dataset_name).(group_name).ICC_mat, 5, 'all'), prctile(ICC_mats.(dataset_name).(group_name).ICC_mat, 95, 'all')];
        % Plot ICC matrix
        subplot(1, 3, i);
        imagesc(ICC_mats.(dataset_name).(group_name).ICC_mat(parc.yeoOrder(ROIs), parc.yeoOrder(ROIs)), range);
        axis square;
        xlabel('ROI test');
        ylabel('ROI retest');
        title(['ICC ', groups_labels{i}]);
        colorbarHandle = colorbar;
        ylabel(colorbarHandle, 'ICC (percentiles 5-95)');
    end
    
    % Add overall title
    if strcmp(dataset_name, 'GENEVA')
        sgtitle(['Figure 3A - ', 'GENEVA']);
    elseif strcmp(dataset_name, 'ADNI')
        sgtitle(['Figure 3B - ', 'ADNI']);
    end
end


% Plot overlap matrix to create Figure 3C
figure;
for i = 1:length(groups)
    group_name = groups{i};
    subplot(1, 3, i);  
    spy(ICC_mats.GENEVA_ADNI.(group_name).overlap(parc.yeoOrder(ROIs), parc.yeoOrder(ROIs)));
    axis square;
    xlabel('ROI test');
    ylabel('ROI retest');
    title(['ICC overlap ', groups_labels{i}]);
    % Calculate number of overlapping edges
    num_overlapping_edge = sum(ICC_mats.GENEVA_ADNI.(group_name).overlap(:));
    subtitle(['n edges= ', num2str(num_overlapping_edge)])
end
sgtitle(['Figure 3C - ', 'GENEVA + ADNI']);

