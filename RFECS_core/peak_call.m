function [pks locs]=peak_call(prob_dist,peak_dist)

[pks locs]=max(prob_dist);

if(length(prob_dist)-locs>peak_dist)
prob_dist2=prob_dist(locs+peak_dist+1:end);
[pks1 locs1]=peak_call(prob_dist2,peak_dist);
pks=[pks pks1];
locs=[locs locs+peak_dist+locs1];
pks1=[];locs1=[];
end

if(locs-1>peak_dist)
prob_dist2=prob_dist(1:locs-peak_dist-1);
[pks1 locs1]=peak_call(prob_dist2,peak_dist);
pks=[pks pks1];
locs=[locs locs1];
pks1=[];locs1=[];
end


end


