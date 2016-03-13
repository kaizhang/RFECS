function forest_predict_par(forest_all,nvar,ntrees,len_var,mod_vec,extn,peak_filt_dist,input_path,output_path,chr_set)
% enhancer predictions using p300 vs TSS-bg

for x=1:length(chr_set)
        chr_set(x)
        all_mods=load(char(strcat(input_path,chr_set(x))));
        all_mods=all_mods(:,mod_vec);
	    index_map=[];
        mean_bg_all=[];
       
        for k=1:len_var
            limit=floor(length(all_mods)/len_var)*len_var+k-1;
            if(limit>length(all_mods))
                limit=limit-len_var;
            end
            index_map=[index_map;[k:len_var:limit]'];
            cd4_set=[];
            for j=1:nvar
                cd4_set=[cd4_set reshape(all_mods(k:limit,j)',len_var,[])'];
            end
            %predicting background score
            y_bg=[];
            parfor i=1:ntrees %for i=1:ntrees
                y_bg(:,i)=str2num(char(eval(forest_all{i},cd4_set)));
            end
            mean_bg=mean(y_bg')';
            mean_bg_all=[mean_bg_all;mean_bg];
        end
        file=strcat(output_path,chr_set(x),extn);
        save(char(file),'mean_bg_all');
        [s sind]=sort(index_map);
        mean_bg_all=mean_bg_all(sind);
        index_sel=find(mean_bg_all>0.5);  
        ind=intersect(find(index_sel>len_var),find(index_sel<=length(index_map)-len_var));
        index_sel=index_sel(ind);
        
        mean_bg_sel=mean_bg_all(index_sel);
        
        patches2=[];
        
        start=index_sel(1);
        prev=index_sel(1);
        for i=2:length(index_sel)
            if(index_sel(i)-prev==1)
                prev=index_sel(i);
            else
                patches2=[patches2;start prev];
                start=index_sel(i);
                prev=start;
            end
        end
        patches2=[patches2;start prev];
        
        patches2_bg={};
        for i=1:length(patches2)
            abc=[];
            for j=patches2(i,1):patches2(i,2)
                abc=[abc mean_bg_all(j)];
            end
            patches2_bg{i}=abc;
        end
        mean_y_bg_uniq=[0:1:ntrees]./ntrees;
        peaks=[];
        peaks_prob=[];
        for i=1:length(patches2)
            vals1=patches2_bg{i};
            [pks locs]=peak_call(vals1,peak_filt_dist); %if peak_dist is 2kb, it should be 20
            peaks=[peaks patches2(i,1)+locs-1]; % change condition to patches2_2kb
            peaks_prob=[peaks_prob pks];
        end
        [s sind]=sort(peaks);
        peaks=peaks(sind);
        peaks_prob=peaks_prob(sind);
        prev=1;
        peaks_filt=[];
        peaks_prob_filt=[];
      
        for i=2:length(peaks)
            if(peaks(i)-peaks(prev)<=peak_filt_dist)
                if(peaks_prob(i)>peaks_prob(prev))
                    prev=i;
                end
            else
                peaks_filt=[peaks_filt;peaks(prev)];
                peaks_prob_filt=[peaks_prob_filt;peaks_prob(prev)];
                prev=i;
            end
        end
        mean_bg_filt=zeros(length(peaks_filt),1);
        for i=1:length(peaks_filt)
            ind1=find(index_sel==peaks_filt(i));
            mean_bg_filt(i)=mean_bg_sel(ind1);
        end
        loc=(peaks_filt+(len_var/2))*100;
filename=strcat(output_path,chr_set(x),extn,'.txt');
fid=fopen(char(filename),'w');
for i=1:length(loc)
fprintf(fid,'%s\t%d\t%f\n',char(chr_set(x)),loc(i),mean_bg_filt(i));
end
end
