%% mdd bd, gmv normative deviation subtyping 
clear
clc

%%% load data
load('hc_dk68_gmv.mat')
load('mdd_bd_dk68_gmv.mat')
load('region_name.mat')
load('age_sex_edu_md_hc.mat')
load('clinical_scores_md.mat') 

disp('loading data finished.....')

%% normative deviation calculation

mean_age1 = mean(age_sex_edu_hc(:,1))
mean_age2 = mean(age_sex_edu_md_data(:,1))
mean_age = (mean_age1+mean_age2)/2

[regional_nd_p,regional_nd_stats] = normative_deviation_hc(age_sex_edu_hc,hc_temp,mean_age);
[regional_nd_md_data,regional_nd_z_md_data] = normative_deviation_md(regional_nd_p,regional_nd_stats,age_sex_edu_md_data,md_temp,mean_age);

%% gmv changes summ data
regional_nd_md_data_1 = zeros(size(md_temp,1),68);  %% z95
regional_nd_md_data_2 = zeros(size(md_temp,1),68);  %% z50
regional_nd_md_data_3 = zeros(size(md_temp,1),68);  %% z05
for i=1:68
    id_temp1 = find(regional_nd_md_data{i}(:,1)>1.96);
    id_temp2 = find(regional_nd_md_data{i}(:,3)<-1.96);
    regional_nd_md_data_1(id_temp1,i) = regional_nd_md_data{i}(id_temp1,1);
    regional_nd_md_data_1(id_temp2,i) = regional_nd_md_data{i}(id_temp2,3);
    regional_nd_md_data_2(:,i) = regional_nd_md_data{i}(:,2);
    regional_nd_md_data_3(:,i) = (regional_nd_md_data{i}(:,1)+regional_nd_md_data{i}(:,2)+regional_nd_md_data{i}(:,3))/3;
end

%% hierarchical clustering
%%% 
Y = zscore(regional_nd_md_data_2');
dis = pdist(Y','euclidean');
dis_squ = squareform(dis);
z = linkage(dis_squ,'ward'); 

%%% color 0.00,0.45,0.74
figure
[tree_p,a,order] = dendrogram(z,0,'ColorThreshold',90);%可视化聚类树

T2=cluster(z,2);%剪枝为三2类

%% subtype 1 and subtype 2,  gmv nd z values 
regional_nd_md_data_1_t1 = regional_nd_md_data_2(T2==1,:);
regional_nd_md_data_1_t2 = regional_nd_md_data_2(T2==2,:);

mean_t1 = mean(regional_nd_md_data_1_t1);
mean_t2 = mean(regional_nd_md_data_1_t2);

mean_all = [mean_t1;mean_t2];

%%% nd 均值偏移；
figure
subplot(2,1,1)
bar(mean_all(:,1:34)')
set(gca,'xtick',[1:34],'xticklabel',regions_name,'XTickLabelRotation',45)
ylabel('Normative Deviation')
subtitle('Left Hemisphere')
subplot(2,1,2)
bar(mean_all(:,35:68)')
set(gca,'xtick',[1:34],'xticklabel',regions_name,'XTickLabelRotation',45)
ylabel('Normative Deviation')
subtitle('Right Hemisphere')

%%% z values
regional_nd_md_data_1_t1 = regional_nd_md_data_1(T2==1,:);
regional_nd_md_data_1_t2 = regional_nd_md_data_1(T2==2,:);

mean_t1 = sum(regional_nd_md_data_1_t1)./sum(sign(abs(regional_nd_md_data_1_t1)));
mean_t2 = sum(regional_nd_md_data_1_t2)./sum(sign(abs(regional_nd_md_data_1_t2)));

mean_t1(isnan(mean_t1)) =0;
mean_t1(isinf(mean_t1)) =0;
mean_t2(isnan(mean_t2)) =0;
mean_t2(isinf(mean_t2)) =0;

sub1_mean_z = mean_t1';
[~,max_id] = max(sub1_mean_z)
regions_name(max_id)

sub2_mean_z = mean_t2';
[~,max_id] = max(sub2_mean_z)
regions_name(max_id-34)


%%%% thres 1
sub1_mean_z_1 = sub1_mean_z;
sub1_mean_z_1(abs(sub1_mean_z_1)<2.3)=0;


sub2_mean_z_1 = sub2_mean_z;
sub2_mean_z_1(abs(sub2_mean_z_1)<2.3)=0;

%%%% thres 2
sub1_mean_z_2 = sub1_mean_z;
sub1_mean_z_2(abs(sub1_mean_z_1)<1.96)=0;


sub2_mean_z_2 = sub2_mean_z;
sub2_mean_z_2(abs(sub2_mean_z_1)<1.96)=0;

%% percentiles of significant gmv nd z scores in all participants of subtype 1 and subtype 2
regional_nd_md_data_1_t1 = regional_nd_md_data_1(T2==1,:);
regional_nd_md_data_1_t2 = regional_nd_md_data_1(T2==2,:);

sign_t1 = sum(regional_nd_md_data_1_t1>0)/sum(T2==1);
sign_t2 = sum(regional_nd_md_data_1_t1<0)/sum(T2==1);

sub1_supra_rate = sign_t1';
sub1_infra_rate = sign_t2';

[~,max_id] = max(sub1_supra_rate)

sign_all = [sign_t1;sign_t2];
sign_all = sign_all*100;

figure
subplot(2,1,1)
bar(sign_all(:,1:34)')
set(gca,'xtick',[1:34],'xticklabel',regions_name,'XTickLabelRotation',45)
ylabel('Subject Number %')
subtitle('Left Hemisphere')
subplot(2,1,2)
bar(sign_all(:,35:68)')
set(gca,'xtick',[1:34],'xticklabel',regions_name,'XTickLabelRotation',45)
ylabel('Subject Number %')
subtitle('Right Hemisphere')
title('Subtype1')
legend('Supra Norm','Infra Norm')

sign_t1 = sum(regional_nd_md_data_1_t2>0)/sum(T2==2);
sign_t2 = sum(regional_nd_md_data_1_t2<0)/sum(T2==2);
sub2_supra_rate = sign_t1';
sub2_infra_rate = sign_t2';

[~,max_id] = max(sub2_infra_rate)

sign_all = [sign_t1;sign_t2];
sign_all = sign_all*100;

figure
subplot(2,1,1)
bar(sign_all(:,1:34)')
set(gca,'xtick',[1:34],'xticklabel',regions_name,'XTickLabelRotation',45)
ylabel('Subject Number %')
subtitle('Left Hemisphere')
subplot(2,1,2)
bar(sign_all(:,35:68)')
set(gca,'xtick',[1:34],'xticklabel',regions_name,'XTickLabelRotation',45)
ylabel('Subject Number %')
subtitle('Right Hemisphere')
legend('Supra Norm','Infra Norm')
title('Subtype2')

%% clinical behaviours 

%%% age differences
scores_md_sub1 = scores_md(T2==1,:);
scores_md_sub2 = scores_md(T2==2,:);

scores_md_age_t1 = scores_md(T2==1,6);
scores_md_age_t2 = scores_md(T2==2,6);

[h,p] = ttest2(scores_md_age_t1,scores_md_age_t2)

%%% duration differences
scores_md_dur_t1 = scores_md(T2==1,11);
scores_md_dur_t2 = scores_md(T2==2,11);
[h,p] = ttest2(scores_md_dur_t1,scores_md_dur_t2)

%%% onset age differences
scores_md_onset_t1 = scores_md_age_t1- scores_md_dur_t1/12;
scores_md_onset_t2 = scores_md_age_t2- scores_md_dur_t2/12;
[h,p] = ttest2(scores_md_onset_t1,scores_md_onset_t2)

%%% clinical symptoms, hamd, hama, bprs, ymrs
scores_md_clin_t1 = scores_md(T2==1,[12,13,14,15]);
scores_md_clin_t2 = scores_md(T2==2,[12,13,14,15]);

for i=1:4
    clin_temp1 = scores_md_clin_t1(:,i);
    clin_temp2 = scores_md_clin_t2(:,i);
    if i>2
        clin_temp1(clin_temp1==0) = [];
        clin_temp2(clin_temp2==0) = [];
    else
        clin_temp1(clin_temp1<7) = [];
        clin_temp2(clin_temp2<7) = [];
    end
    [h(i),p(i)] = ttest2(clin_temp1,clin_temp2);
end

%%% diagnosis mdd, bd
scores_md_diag_t1 = scores_md(T2==1,[2]);
scores_md_diag_t2 = scores_md(T2==2,[2]);

disp('-------------mdd------------')
sum(scores_md_diag_t1(:,1)==2)
sum(scores_md_diag_t2(:,1)==2)
disp('-------------bd------------')
sum(scores_md_diag_t1(:,1)==4)
sum(scores_md_diag_t2(:,1)==4)

scores_md_diag_subtype=[];
scores_md_diag_subtype(1,1) = sum(scores_md_diag_t1(:,1)==2);
scores_md_diag_subtype(1,2) = sum(scores_md_diag_t2(:,1)==2);
scores_md_diag_subtype(2,1) = sum(scores_md_diag_t1(:,1)==4);
scores_md_diag_subtype(2,2) = sum(scores_md_diag_t2(:,1)==4);
bar(scores_md_diag_subtype,'BarLayout','stacked')
set(gca,'XTicklabel',{'Subtype1','Subtype2'});
legend('MDD','BD')














