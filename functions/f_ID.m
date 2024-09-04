function [idiff_m, idiff_se, idiff_vect,  success, iself_m, iself_sd,Iself_vect, Iothers_id,Iothers_m ] = f_ID(ID)

%ID matrix and indexes

N_subj = size(ID,1);
mask_diag = logical(eye(N_subj));
%I-Self
Iself_vect = (ID(mask_diag));
%Calculate Iothers mean for each subject
Iothers_id = zeros(N_subj,1);
for i = 1:N_subj
    Iothers_id(i)=(sum(ID(i,:))+ sum(ID(:,i))-2*ID(i,i))/(2*N_subj-2);
end
%Calculate success rate and i diff
success = (nnz(Iself_vect > Iothers_id))*100/N_subj;
idiff_vect = Iself_vect - Iothers_id;
idiff_m = mean(Iself_vect - Iothers_id);
idiff_sd = std(Iself_vect - Iothers_id);
idiff_se = std(Iself_vect - Iothers_id)/sqrt(length(Iself_vect - Iothers_id));

iself_m = mean(Iself_vect);
iself_sd = std(Iself_vect);
Iothers_m = mean(Iothers_id);

