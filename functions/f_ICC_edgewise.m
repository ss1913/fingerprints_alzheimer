function [ICC_vec, ICC_mat] = f_ICC_edgewise(Atest,Aretest,N)
% Evaluate ICC on the connICA matrices, (Parfor to Speed up computation) 
% Add ICC & Network Scripts to path
addpath(fullfile(pwd,'..','functions/ICC/'))

%% First plot: mean triu ICC vs num ICA comps
numEdges = size(Atest,1);

ICC_vec = zeros(1,numEdges);        
ICC_mat = zeros(N,N);
mask_ut = triu(true(N),1);
%parpool;
for comp=1:numEdges
    %data4icc_A = data4icc(:,:,comp);
    data4icc_A = [Atest(comp,:)' Aretest(comp,:)'];
    rows2delete_A = isnan(sum(data4icc_A,2));
    data4icc_A(rows2delete_A,:) = [];
    ICC_vec(1,comp) = ICC(data4icc_A,'1-1') ; 
end
ICC_mat(mask_ut) = ICC_vec;
ICC_mat = ICC_mat + ICC_mat';
%delete(gcp);
return;