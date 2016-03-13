%Calculate OOB variable importance for forest
function [misclass_perm var_imp]=OOBVarImp(forest_p300_tssbg,nvar,rand_set_all,set_p300_tssbg,vals,len_vars)
% forest_p300_tssbg is the trained forest,nvar is number of histone marks,rand_set_all is the matrix containing index of elements chosen in each forest,set_p300_tssbg is the training set , vals in the labels of elements in training set,len_vars is the vector containing length of each element in the training set

%misclassification probabilities for each variables in each tree
ntrees=length(forest_p300_tssbg);
cum_len=cumsum(len_vars);

for i=1:ntrees
    oob_index=setdiff([1:length(rand_set_all)],rand_set_all(:,i));
    ts3=[];y_oob=[];
    %ts3=zeros(length(oob_index)*26,nvar*20);
    %ts3(1:length(oob_index),:)=set_p300_tssbg(oob_index,:);
    ts3=set_p300_tssbg(oob_index,:);
    
    for j=1:nvar
        rand_col=randperm(length(oob_index));
        
        if(j>1)
        ind_len=[cum_len(j-1)+1:cum_len(j)];
        var_col=set_p300_tssbg(oob_index,ind_len);
        var_col2=var_col(rand_col,:);            
        ts3=[ts3;set_p300_tssbg(oob_index,1:cum_len(j-1)) var_col2 set_p300_tssbg(oob_index,cum_len(j)+1:end)];
        else
        ind_len=[1:cum_len(j)];
        var_col=set_p300_tssbg(oob_index,ind_len);
        var_col2=var_col(rand_col,:);            
        ts3=[ts3;var_col2 set_p300_tssbg(oob_index,cum_len(j)+1:end)];
        end
    end
    y_oob=str2num(char(eval(forest_p300_tssbg{i},ts3)));
    
    for j=1:nvar+1
        misclass_perm(i,j)=length(find(round(y_oob(length(oob_index)*(j-1)+1:length(oob_index)*j))~=        vals(oob_index)))/length(oob_index);
    end
    
end

% find difference in misclasification errors

diff_var=misclass_perm(:,2:end)-repmat(misclass_perm(:,1),1,nvar);

mean_diff=mean(diff_var);
std_diff=std(diff_var);

var_imp=mean_diff./std_diff;

