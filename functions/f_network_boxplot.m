function [Y_within, Y_between, chisq_table ] = f_network_boxplot(ICC_mat,Group, yeoROIs,key_WithinNetworks,key_BetweenNetworks)

N = size(ICC_mat,1); 
n_yeo = 8;
yeoROIs = yeoROIs(1:N,1);

count = 1;
for i=1:n_yeo
    for j=1:n_yeo
        if i == j
            aux{i,1}= ICC_mat(yeoROIs==i,yeoROIs==j); %within networks 
        else  
            aux2{i,:} = ICC_mat(yeoROIs==i,yeoROIs(:)); %between networks
        end
        count = count + 1;
    end
end


%% Compute average of within-networks edges

aux_reshape = [];
grp_reshape = [];
for i =1:size(aux,1)
    a = reshape(aux{i,1},[],1)';
    aux_reshape  = [aux_reshape, a];
    b = i*ones(1,length(reshape(aux{i,1},[],1)));
    grp_reshape  = [grp_reshape, b];
end

T = key_WithinNetworks;
for j =1:size(T,1)
    for i = 1:size(grp_reshape,2)
        if grp_reshape(i)== T.Var1(j)
            New_grp_reshape(i,1)=T.Var2(j);
            class_reshape(i,2)=T.Var3(j);
        end
    end

end
grp_reshape = New_grp_reshape;
class_reshape = class_reshape(:,2);
gr = strings(1, size(grp_reshape,1));
gr(:) = Group;
type = strings(1, size(grp_reshape,1));
type(:)= "Within";
data = table(aux_reshape', grp_reshape,class_reshape,gr',type','VariableNames', {'ICC', 'Network', 'Class', 'Group', 'Type'});
statarray = grpstats(data,"Network", "mean", DataVars="ICC");
statarray_se = grpstats(data,"Network", "sem", DataVars="ICC");
Y_within = statarray.mean_ICC*100;
Y_within_err= statarray_se.sem_ICC*100;
X_within = categorical(statarray.Network);

%% Compute average of between networks edges

aux2_reshape = [];
grp2_reshape = [];

for i =1:size(aux2,1)
    a = reshape(aux2{i,1},[],1)';
    aux2_reshape  = [aux2_reshape, a];
    b = i*ones(1,length(reshape(aux2{i,1},[],1)));
    grp2_reshape  = [grp2_reshape, b];
end
T = key_BetweenNetworks;
for j =1:size(T,1)
    for i = 1:size(grp2_reshape,2)
        if grp2_reshape(i)== T.Var1(j)
            New_grp2_reshape(i,1)=T.Var2(j);
            class2_reshape(i,2)=T.Var3(j);
        end
    end
 end

grp2_reshape = New_grp2_reshape;
class2_reshape = class2_reshape(:,2);
gr2 = strings(1, size(grp2_reshape,1));
gr2(:) = Group;
type2 = strings(1, size(grp2_reshape,1));
type2(:)= "Between";
data2 = table(aux2_reshape', grp2_reshape,class2_reshape,gr2',type2','VariableNames', {'ICC', 'Network', 'Class', 'Group', 'Type'});
statarray2 = grpstats(data2,"Network", "mean", DataVars="ICC");
statarray2_se = grpstats(data2,"Network", "sem", DataVars="ICC");
Y_between = statarray2.mean_ICC*100;
Y_between_err= statarray2_se.sem_ICC*100;
X_between = categorical(statarray2.Network);

%% Combine within and between networks and make graph

for ii = numel(X_within):-1:1
    X_wb(2*ii-[0 1]) = [X_between(ii),X_within(ii)];
end

for ii = numel(Y_within):-1:1
    Y_wb(2*ii-[0 1]) = [Y_between(ii),Y_within(ii)];
end

for ii = numel(Y_within_err):-1:1
    Y_wb_err(2*ii-[0 1]) = [Y_between_err(ii),Y_within_err(ii)];
end
X_wb = reordercats(X_wb,{'VIS','VIS-between','SMT','SMT-between','DA','DA-between','SA','SA-between','L','L-between','FPN','FPN-between','DMN','DMN-between','SBC','SBC-between'});

figure;
bar(X_wb, Y_wb); title(Group);ylim([0 35])
hold on
errorbar(X_wb, Y_wb,Y_wb_err, '.' )


%% Compare within and between-networks number of edges using Chi-Square test 

% get unique groups

group_vars = unique(data.Network);
group_vars2 = unique(data2.Network);

chisq_table = table(group_vars, 'VariableNames', {'Networks'});
chisq_table.p_value = zeros(size(group_vars));
chisq_table.chi_statistic = zeros(size(group_vars));

% iterate over groups and perform Chi-Square test
for i = 1:numel(group_vars)
    % Get the data for the contingency table
    ICC1 = data.ICC(strcmp(data.Network, group_vars{i}))'; %vector of ICC within network 
    ICC2 = data2.ICC(strcmp(data2.Network, group_vars2{i}))'; %vector of ICC between network
    cond1 = zeros(1,(size(ICC1,2)));
    cond2 = ones(1,(size(ICC2,2)));
    ICC= [ICC1, ICC2];
    cond = [cond1, cond2];

    % Perform the Chi-square test
    [tbl,chi2stat,pval] = crosstab(ICC,cond);
    
    % Store the results in the ttest_table
    chisq_table .p_value(i) = pval;
    chisq_table.chi_statistic (i) = chi2stat;

end

% Reorder the rows of ttest_table according to the desired sequence
Network_order = {'VIS','SMT','DA','SA','L','FPN','DMN','SBC'}; 
[~, idx] = ismember(Network_order, chisq_table.Networks);
chisq_table = chisq_table(idx, :);









