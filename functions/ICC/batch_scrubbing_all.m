%% Batch to load FC according to the table

Path2data = fullfile('..','data_full'); 
Path2Table = fullfile('..','tables', 'all/'); 
Path2Results = fullfile('..','results', 'all/'); 
T = readtable(fullfile(Path2Table,'all_20210126.csv'));
%% Read ID from Table, and select FC folder accordingly
N_subj = size(T,1);
N = 278;

Scrub_vect = zeros(1,N_subj);
Scrub_vect_test = zeros(1,N_subj);
Scrub_vect_retest = zeros(1,N_subj);
for s=1:N_subj
    disp('processing subject');
    Subj_ID = string(T.short_core_id(s));
    %Subj_ID = num2str(T.short_core_id(s));
    disp(Subj_ID);
    fmri_folder = T.T_0_fMRI_folder{s};
    Path2scrub = fullfile(Path2data,Subj_ID,fmri_folder,'GSreg_yes','6_scrubbing.mat');
    load(Path2scrub) % 
    Scrub_vect(s) = nnz(scrubbing(1:200)==0); %scrubbing zero means the volume should be scrubbed.
    Scrub_vect_test(s) = nnz(scrubbing(1:100)==0);
    Scrub_vect_retest(s) = nnz(scrubbing(101:200)==0);
end


save(fullfile(Path2Results,'Scrubbing_20210126_all.mat'),'Scrub_vect','Scrub_vect_test','Scrub_vect_retest');