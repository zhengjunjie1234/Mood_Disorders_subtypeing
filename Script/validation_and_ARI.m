%%
clear
clc

load('discovery_dataset.mat')
%% ARI test

Y = zscore(regional_nd_md_data_2');
dis = pdist(Y','euclidean');
dis_squ = squareform(dis);
z = linkage(dis_squ,'ward'); 

%%% color 0.00,0.45,0.74
figure
[tree_p,a,order] = dendrogram(z,0,'ColorThreshold',90);%可视化聚类树

T2=cluster(z,2);%剪枝为三2类

Kmax = 6;
clust = zeros(size(Y,2),Kmax);
for k=1:Kmax
    clust(:,k) = cluster(z, 'maxclust', k);
end

%%
%%% 10 fold 
fold_10_sample ={};
fold_10_sample_1 = {};
for i=1:5000
    temp_id = randsample([1:174],18);
    fold_10_sample{i} = setdiff([1:174],temp_id);
    fold_10_sample_1{i} = clust(fold_10_sample{i},:);
end

ari = []
for i=1:5000
    Y_temp = Y(:,fold_10_sample{i});
    dis = pdist(Y_temp','euclidean');
    dis_squ = squareform(dis);
    z = linkage(dis_squ,'ward'); 
    clust_temp = zeros(size(Y_temp,2),Kmax);
    for k=1:Kmax
        clust_temp(:,k) = cluster(z, 'maxclust', k);
    end
    
    for k=2:Kmax
        ari(i,k) = pairwiseindex(clust_temp(:,k),fold_10_sample_1{i}(:,k));
    end
end

ari_mean = mean(ari)
ari_std = std(ari)
%% plot




%%
