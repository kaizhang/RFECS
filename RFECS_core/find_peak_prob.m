function prob=find_peak_prob(values,prob_mat)
mid=floor(1+length(values))/2;
prob_mid=floor((1+size(prob_mat,2))/2);
prob=prob_mat(values(mid),prob_mid);
for i=1:length(values)
    if(i~=mid)
prob=prob*prob_mat(values(i),-mid+prob_mid+i);
    end
end

