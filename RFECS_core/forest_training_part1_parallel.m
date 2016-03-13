function forest_training_part1_parallel(set_p300_tssbg,vals,nvar,ntrees,len_var,forest_file)
parfor i=1:ntrees
rand_set_all(:,i)=randi(length(vals),[length(vals),1]);
forest_p300_tssbg_all{i}=multiclasstree(set_p300_tssbg(rand_set_all(:,i),:),vals(rand_set_all(:,i)),'categorical',[1:len_var*nvar],'categorical_group',nvar,'nvartosample',ceil(sqrt(nvar)));
end

save(char(forest_file),'forest_p300_tssbg_all','rand_set_all');
