%% Batch for Supplementary 
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

% Load outputs
load(fullfile(Path2Outputs, 'ICC_mats.mat'));

%% Define constants and variables
% Constants
groups = {'CU_minus', 'MCI_plus', 'Dementia_plus'};
datasets = {'GENEVA', 'ADNI'};
N = size(parc.yeoROIs, 1);


% Set FSL environment
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
path2FSL = '/usr/local/fsl/bin'; % change this according to your FSL directory

%% Theshold maps and binarize them 

for j = 1:length(datasets)
    dataset_name = datasets{j};
    for i = 1:length(groups)
        group_name = groups{i};

        % Threshold map
        thr = ICC_mats.(dataset_name).(group_name).nodalstrenght_percentiles_75_100(1,1);
        fileIn = fullfile(Path2Outputs, sprintf("brain_renders/ICC_sig_%s_%s_shen_dil.nii.gz", dataset_name, group_name));
        fileOut = fullfile(Path2Outputs, sprintf("brain_renders/ICC_sig_%s_%s_shen_dil_bin.nii.gz", dataset_name, group_name));
        sentence = sprintf('%s/fslmaths %s -thr %s %s',path2FSL,fileIn, thr, fileOut);
        [status,result] = system(sentence); 
        
        % Binarize map
        sentence = sprintf('%s/fslmaths %s -bin %s',path2FSL,fileOut,fileOut);
        [status,result] = system(sentence);
    end
end

%% Generate a trasnfromation matrix to bring nifti to MNI template

% Create 'neurosynth' folder in outputs 
mkdir(fullfile(Path2Outputs,"neurosynth"));

fileIn = fullfile(Path2Inputs,"shen_dil.nii.gz");
fileRef = fullfile(Path2Inputs,"Shen_MNI_2mm.nii.gz");
fileOut = fullfile(Path2Outputs,"neurosynth/shen2MNI.mat");

sentence = sprintf('%s/flirt -in %s -ref %s -omat %s ', path2FSL, fileIn, fileRef, fileOut);
[status, result] = system(sentence);

%% Combine two cohorts 

for i = 1:length(groups)
        group_name = groups{i};

        % Combine the two cohorts maps
        fileIn1 = fullfile(Path2Outputs, sprintf("brain_renders/ICC_sig_GENEVA_%s_shen_dil_bin.nii.gz", group_name));
        fileIn2 = fullfile(Path2Outputs, sprintf("brain_renders/ICC_sig_ADNI_%s_shen_dil_bin.nii.gz", group_name));
        fileOut = fullfile(Path2Outputs, sprintf("neurosynth/ICC_sig_GENEVA_ADNI_%s_shen_dil_bin.nii.gz", group_name));
        sentence = sprintf('%s/fslmaths %s -add %s %s',path2FSL,fileIn1, fileIn2, fileOut);
        [status,result] = system(sentence);
        
        % Binarize maps
        sentence = sprintf('%s/fslmaths %s -thr 2 %s',path2FSL,fileOut, fileOut);
        [status,result] = system(sentence);

        % Register onto MNI2mm
        fileIn = fullfile(Path2Outputs, sprintf("neurosynth/ICC_sig_GENEVA_ADNI_%s_shen_dil_bin.nii.gz", group_name));
        fileRef = fullfile(Path2Inputs,"Shen_MNI_2mm.nii.gz");
        fileMat = fullfile(Path2Outputs,"/neurosynth/shen2MNI.mat");       
        fileOut = fullfile(Path2Outputs, sprintf("neurosynth/ICC_sig_GENEVA_ADNI_%s_shen_dil_bin.nii.gz", group_name));
        sentence = sprintf('%s/flirt -in %s -ref %s -applyxfm -init %s -out %s',path2FSL,fileIn, fileRef, fileMat, fileOut);
        [status,result] = system(sentence);

        % Threshold maps
        fileIn = fileOut;
        fileOut = fullfile(Path2Outputs, sprintf("neurosynth/ICC_sig_GENEVA_ADNI_%s_shen_dil_bin.nii.gz", group_name));
        sentence = sprintf('%s/fslmaths %s -thr 1 %s',path2FSL,fileIn, fileOut);
        [status,result] = system(sentence);

        % Binarize maps
        fileIn = fileOut;
        fileOut = fullfile(Path2Outputs, sprintf("neurosynth/ICC_sig_GENEVA_ADNI_%s_shen_dil_bin.nii.gz", group_name));
        sentence = sprintf('%s/fslmaths %s -bin %s',path2FSL,fileIn, fileOut);
        [status,result] = system(sentence);

end












