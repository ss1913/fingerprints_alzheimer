function [idiff_null, success_null] = f_ID_permutation(ID,null_runs)

%Permutation testing Identifiability matrix
N_subj = size(ID,1);
mask_diag = logical(eye(N_subj));
success_null = zeros(1,null_runs);
idiff_null = zeros(1,null_runs);

for run = 1:null_runs
    
    ID_perm = ID(randperm(N_subj),randperm(N_subj));
    %I-Self 
    Iself_vect = (ID_perm(mask_diag));
    %Calculate Iothers mean for each subject
    Iothers_id = zeros(N_subj,1);
    for i = 1:N_subj
        Iothers_id(i)=(sum(ID_perm(i,:))+ sum(ID_perm(:,i))-2*ID_perm(i,i))/(2*N_subj-2);
    end 
    %Calculate success rate and i diff
    success_null(run) = (nnz(Iself_vect > Iothers_id))*100/N_subj;
    idiff_null(run) = mean(Iself_vect - Iothers_id);


end
