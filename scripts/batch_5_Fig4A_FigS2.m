%% Batch for Supplementary Figure 2 and Figure 4A
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

% Load Data
load(fullfile(Path2Outputs, 'ICC_mats.mat'));

% Load parcellation
parc = load(fullfile(Path2Inputs, 'shen_yeo_RS7.mat'));

%% Define constants and variables

% Constants
N = size(parc.yeoROIs, 1);
groups = {'CU_minus', 'MCI_plus', 'Dementia_plus'};
groups_labels = {'CU AB-', 'MCI AB+', 'AD Dementia'};
datasets = {'GENEVA', 'ADNI'};

%% Supplementary Figure 2

% Calculate average of edges in overlap ICC matric per each within and between network to create Supplementary Figure 2 
for i = 1:length(groups)
    group_name = groups{i};
    [within, between, chisq_table] = f_network_boxplot(ICC_mats.GENEVA_ADNI.(group_name).overlap, group_name, parc.yeoROIs, parc.key_WithinNetworks, parc.key_BetweenNetworks);
    ICC_mats.GENEVA_ADNI.(group_name).within = within;
    ICC_mats.GENEVA_ADNI.(group_name).between = between;
    ICC_mats.GENEVA_ADNI.(group_name).chisquaretable= chisq_table;
end 

%% Figure 4A

% Compute ratio within/between with respect to CU_minus to create Figure 4A
metrics = {'within', 'between'};
label = categorical({'VIS','SMT','DA','SA','L','FPN','DMN','SBC'}); 

% Loop over each group and each metric
for j = 1:length(metrics)
    metric_name = metrics{j};
    
    for i = 2:length(groups)
        group_name = groups{i};
        
        % Determine the healthy reference (CU_minus) for the current metric
        healthy_ref = ICC_mats.GENEVA_ADNI.CU_minus.(metric_name);
        
        % Calculate the ratio disease/health
        current = ICC_mats.GENEVA_ADNI.(group_name).(metric_name);
        ratio = (current./healthy_ref) - 1;
        ICC_mats.GENEVA_ADNI.(group_name).(sprintf('ratio_%s', metric_name))= ratio;
    end

end

mean(ICC_mats.GENEVA_ADNI.MCI_plus.ratio_within)
mean(ICC_mats.GENEVA_ADNI.MCI_plus.ratio_between)
mean(ICC_mats.GENEVA_ADNI.Dementia_plus.ratio_within)
mean(ICC_mats.GENEVA_ADNI.Dementia_plus.ratio_between)

% Combine data for figure 4
tbl = table(label', ...
    ICC_mats.GENEVA_ADNI.Dementia_plus.ratio_between,...
    ICC_mats.GENEVA_ADNI.Dementia_plus.ratio_within, ...
    ICC_mats.GENEVA_ADNI.MCI_plus.ratio_between,...
    ICC_mats.GENEVA_ADNI.MCI_plus.ratio_within); 

% Reorder networks
X = reordercats(categorical(tbl{1:8,1}), {'SBC','DMN','FPN','L','SA','DA','SMT','VIS'});

figure; 
barh(X, tbl{1:8,2:5}, 'BarWidth', 0.8);

% Manually specify legend labels
legendLabels = {'Dementia Plus - Between','Dementia Plus - Within'...
    'MCI Plus - Between','MCI Plus - Within'};
hLegend = legend(legendLabels, 'Location', 'best');

xlim([-1, 12])

sgtitle(['Figure 4A - ', 'Ratio Disease/Health']);


%% Create a table with chisquare test within-between networks for Figure 4B

for i = 1:length(groups)
    group_name = groups{i};
    % Apply Bonferroni correction
    ICC_mats.GENEVA_ADNI.(group_name).chisquaretable.p_value_bonfcorr = ICC_mats.GENEVA_ADNI.(group_name).chisquaretable.p_value.*size(ICC_mats.GENEVA_ADNI.(group_name).chisquaretable,1);       
end

chisq_table = table(ICC_mats.GENEVA_ADNI.CU_minus.chisquaretable.Networks, ...
    ICC_mats.GENEVA_ADNI.CU_minus.chisquaretable.chi_statistic, ...
    ICC_mats.GENEVA_ADNI.MCI_plus.chisquaretable.chi_statistic, ...
    ICC_mats.GENEVA_ADNI.Dementia_plus.chisquaretable.chi_statistic, ...
    'VariableNames', {'Network', 'CU_minus', 'MCI_plus', 'Dementia_plus'});

writetable(chisq_table, fullfile(Path2Outputs,'chi_square_within_between.csv'));

% For Figure 4B see batch_6_Fig_4B.R
