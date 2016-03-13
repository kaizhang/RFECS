% first part of training to generate forest for classification
function forest_training_part1(set_p300_tssbg,vals,nvar,ntrees,len_var,forest_file)
% set_p300_tssbg is the training set composed of p300/TSS/random background in the ratio 1:1:x, vals is the labels of the training set with ones corresponding to p300 and zeros to the rest,nvar is number of histone modifications,len_var is length of histone modification vector,ntrees is number of trees to be grown,forest_file is string containing name of mat-file where prob_bg_p300_dist should be saved
for i=1:ntrees
rand_set_all(:,i)=randint(1,length(vals),[1 length(vals)]);
forest_p300_tssbg_all{i}=multiclasstree(set_p300_tssbg(rand_set_all(:,i),:),vals(rand_set_all(:,i)),'categorical',[1:len_var*nvar],'categorical_group',nvar,'nvartosample',ceil(sqrt(nvar)));
end

save(char(forest_file),'forest_p300_tssbg_all','rand_set_all');
