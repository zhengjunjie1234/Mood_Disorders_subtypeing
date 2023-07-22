%%
clear
clc

%% load data
[clinical_data,head] = xlsread('G:\shenyang_methylation\mdd_bd_gmscn_subtype\clinic_info.xlsx');

hc_id = xlsread('G:\shenyang_methylation\mdd_bd_gmscn_subtype\omics\hc.xlsx');
subtype1_id = xlsread('G:\shenyang_methylation\mdd_bd_gmscn_subtype\omics\subtype1.xlsx');
subtype2_id = xlsread('G:\shenyang_methylation\mdd_bd_gmscn_subtype\omics\subtype2.xlsx');

%% education
edu = xlsread('G:\shenyang_methylation\mdd_bd_gmscn_subtype\education.xlsx');
edu(isnan(edu)) = 0;

clinical_data(isnan(clinical_data))=0;
clinical_data(clinical_data(:,3)>10,3)=4;

[~,id] = ismember(subtype1_id,clinical_data(:,1));
subtype1_clin = clinical_data(id,:);
subtype1_edu = edu(id,2);

[~,id] = ismember(subtype2_id,clinical_data(:,1));
subtype2_clin = clinical_data(id,:);
subtype2_edu = edu(id,2);

[~,id] = ismember(hc_id,clinical_data(:,1));
hc_clin = clinical_data(id,:);
hc_edu = edu(id,2);

%% wsct scores regress out education 

subtype1_clin_wsct = subtype1_clin(:,39:43);
subtype2_clin_wsct = subtype2_clin(:,39:43);
hc_wsct = hc_clin(:,39:43); 

subtype1_edu1 = subtype1_edu;
subtype1_edu1(subtype1_clin_wsct(:,1)==0) = [];

subtype2_edu1 = subtype2_edu;
subtype2_edu1(subtype2_clin_wsct(:,1)==0) = [];

hc_edu1 = hc_edu;
hc_edu1(hc_wsct(:,1)==0) = [];

subtype1_clin_wsct(subtype1_clin_wsct(:,1)==0,:) = [];
subtype2_clin_wsct(subtype2_clin_wsct(:,1)==0,:) = [];
hc_wsct(hc_wsct(:,1)==0,:) = [];

clin_wsct_u = [subtype1_clin_wsct;subtype2_clin_wsct;hc_wsct];
edu_u = [subtype1_edu1;subtype2_edu1;hc_edu1];

for i=1:5
    [bb,dev,stats] = glmfit(edu_u,clin_wsct_u(:,i));
    clin_wsct_u_r(:,i)= clin_wsct_u(:,i) - edu_u*bb(2);
end

subtype1_clin_wsct1 = clin_wsct_u(1:length(subtype1_edu1),:);
subtype2_clin_wsct1 = clin_wsct_u(length(subtype1_edu1)+1:length(subtype1_edu1)+length(subtype2_edu1),:);
hc_clin_wsct1 = clin_wsct_u(length(subtype1_edu1)+length(subtype2_edu1)+1:end,:);

%%%% statistical and plot
mean(hc_clin_wsct1)
mean(subtype1_clin_wsct1)
mean(subtype2_clin_wsct1)

std(hc_clin_wsct1)
std(subtype1_clin_wsct1)
std(subtype2_clin_wsct1)

h=[];p=[],t=[];
[h,p(1,:),c,s] = ttest2(subtype1_clin_wsct1,hc_clin_wsct1);
t(1,:) = s.tstat;
[h,p(2,:),c,s] = ttest2(subtype2_clin_wsct1,hc_clin_wsct1);
t(2,:) = s.tstat;
[h,p(3,:),c,s] = ttest2(subtype1_clin_wsct1,subtype2_clin_wsct1);
t(3,:) = s.tstat;

p_fdr = mafdr(p(:),'BHFDR',true);
p_fdr = reshape(p_fdr,3,5);

%%% plot
mean_values = [];std_values=[];
mean_values(1,:) = mean(subtype1_clin_wsct1);
mean_values(2,:) = mean(subtype2_clin_wsct1);
mean_values(3,:) = mean(hc_clin_wsct1);

std_values(1,:) = std(subtype1_clin_wsct1)/5;
std_values(2,:) = std(subtype2_clin_wsct1)/5;
std_values(3,:) = std(hc_clin_wsct1)/5;
std_values(4,:) = zeros([1,5]);

wsct_name = head(39:43)
wsct_name_new={};
for i=1:5
    wsct_name_new{i} = wsct_name{i}(6:end);
end

%%% plot
h = figure('position',[200,200,800,600])
a = bar(mean_values')
a(1).FaceColor = [0.85,0.33,0.10];
a(2).FaceColor = [0.00,0.45,0.74];

hold on
errorbar([1:5]-0.225,mean_values(1,:),std_values(1,:) ,'k' , 'Linestyle', 'None','Linewidth',2);
hold on
errorbar([1:5],mean_values(2,:),std_values(2,:), 'k' , 'Linestyle', 'None','Linewidth',2);
hold on
errorbar([1:5]+0.225,mean_values(3,:),std_values(3,:), 'k' , 'Linestyle', 'None','Linewidth',2);
hold on
scatter([1:4]-0.225,mean_values(1,1:4)+std_values(1,1:4).*[2,4,2,2],'rx','Linewidth',2)
set(gca,'xtick',[1:5],'xticklabel',wsct_name_new,'XTickLabelRotation',45,'FontWeight','bold','FontSize',14)
legend('Subtype 1','Subtype 2','hc')
ylabel('WCST Scores','FontWeight','bold','FontSize',14)


%% hamd items
[~,id] = ismember(subtype1_id,hamd(:,1));
subtype1_clin_hamd_t = hamd(id,6);

subtype1_clin_hamd = subtype1_clin(:,2:18);
subtype1_clin_hamd_f = subtype1_clin(:,19:22);
subtype1_clin_hamd_s = sum(subtype1_clin_hamd_f,2);

[~,id] = ismember(subtype2_id,hamd(:,1));
subtype2_clin_hamd_t = hamd(id,6);

subtype2_clin_hamd = subtype2_clin(:,2:18);
subtype2_clin_hamd_f = subtype2_clin(:,19:22);
subtype2_clin_hamd_s = sum(subtype2_clin_hamd_f,2);

subtype1_clin_hamd(subtype1_clin_hamd_t<7,:) = [];
subtype2_clin_hamd(subtype2_clin_hamd_t<7,:) = [];
subtype1_clin_hamd_t(subtype1_clin_hamd_t<7) = [];
subtype2_clin_hamd_t(subtype2_clin_hamd_t<7) = [];


[h,p] = ttest2(subtype1_clin_hamd_t,subtype2_clin_hamd_t)

[h,p] = ttest2(subtype1_clin_hamd,subtype2_clin_hamd)
hamd_factor_name = head(2:18);
hamd_factor_name(p<0.05)

mean_values = [mean(subtype1_clin_hamd);mean(subtype2_clin_hamd)];
std_values = [std(subtype1_clin_hamd);std(subtype2_clin_hamd)];
std_values = std_values/5;
std_values(3,:) = zeros([1,17]);
h = figure('position',[200,200,800,600])
a = bar(mean_values')
a(1).FaceColor = [0.85,0.33,0.10];
a(2).FaceColor = [0.00,0.45,0.74];
hold on
errorbar([1:17]-0.15,mean_values(1,:),std_values(3,:),std_values(1,:) ,'k' , 'Linestyle', 'None','Linewidth',1);
hold on
errorbar([1:17]+0.15,mean_values(2,:),std_values(3,:),std_values(2,:), 'k' , 'Linestyle', 'None','Linewidth',1);
hold on
scatter(find(p<0.05),mean_values(2,find(p<0.05))+1,'rx','Linewidth',2)
set(gca,'xtick',[1:17],'xticklabel',hamd_factor_name,'XTickLabelRotation',45,'FontWeight','bold','FontSize',14)
ylabel('HAMD scores','FontWeight','bold','FontSize',14)
legend('Subtype 1','Subtype 2')

%% clinical symptoms network (hamd items , hama items, bprs factors )
head_new = head;
head_id = setdiff([1:58],[1,19:22,37:53]);
head_new([1,19:22,37:53]) = [];

a = sum(subtype1_clin(:,2:18),2);
over_sub1 = subtype1_clin(find(a>=7),head_id);
a = sum(subtype2_clin(:,2:18),2);
over_sub2 = subtype2_clin(find(a>=7),head_id);

sub1_clin_net = zeros(length(head_new),length(head_new));

for i=1:1000
    rand_s = randsample([1:size(over_sub1,1)],size(over_sub2,1));
    over_sub_r = over_sub1(rand_s,:);
    sub1_clin_net = sub1_clin_net + corr(over_sub1);
end
sub1_clin_net = sub1_clin_net/1000;
sub2_clin_net = corr(over_sub2);

sub1_clin_net = abs(sub1_clin_net);
sub2_clin_net = abs(sub2_clin_net);

sub1_clin_net_thr = sub1_clin_net;
sub2_clin_net_thr = sub2_clin_net;

sub1_clin_net_thr(sub1_clin_net_thr<0.2) = 0;
sub2_clin_net_thr(sub2_clin_net_thr<0.2) = 0;

%% clinical symptoms network properties
%%% degree
node_s_sub1 = sum(sign(sub1_clin_net_thr))-1;
node_s_sub2 = sum(sign(sub2_clin_net_thr))-1;
node_s_de = [node_s_sub1;node_s_sub2]';

%%% strength
node_s_sub1 = sum((sub1_clin_net_thr))-1;
node_s_sub2 = sum((sub2_clin_net_thr))-1;
node_s_str = [node_s_sub1;node_s_sub2]';

node_s_sub1 = betweenness_wei(sub1_clin_net_thr);
node_s_sub2 = betweenness_wei(sub2_clin_net_thr);
node_s_bc = [node_s_sub1,node_s_sub2];

node_name = head_new

for i=1:17
    node_name{i} = strcat(node_name{i},'(HAMD)');
end

for i=18:31
    if i~=23
        node_name{i} = strcat(node_name{i},'(HAMA)');
    end
end

sub1_clin_net_thr_d = sub1_clin_net;
sub2_clin_net_thr_d = sub2_clin_net;
sub1_clin_net_thr_d(sub1_clin_net_thr_d<0.2) = 0;
sub2_clin_net_thr_d(sub2_clin_net_thr_d<0.2) = 0;

d1 = density_und(sub1_clin_net_thr_d)  %%% BCT toolkit
d2 = density_und(sub2_clin_net_thr_d)

s1 = mean(strengths_und(sub1_clin_net_thr_d))  %%% BCT toolkit
s2 = mean(strengths_und(sub2_clin_net_thr_d))


density_thr = [0.1:0.01:0.2];
sub1_clin_net_thr_d = sub1_clin_net;
sub2_clin_net_thr_d = sub2_clin_net;
density_sub1 = [];
density_sub2 = [];
for i=1:length(density_thr)
    sub1_clin_net_thr_d = sub1_clin_net;
    sub2_clin_net_thr_d = sub2_clin_net;
    sub1_clin_net_thr_d(sub1_clin_net_thr_d<density_thr(i)) = 0;
    sub2_clin_net_thr_d(sub2_clin_net_thr_d<density_thr(i)) = 0;
    density_sub1(i) = density_und(sub1_clin_net_thr_d);
    density_sub2(i) = density_und(sub2_clin_net_thr_d);
end

%%% strength 
density_thr_1 = [0.1:0.05:0.6];
sub1_clin_net_thr_d = sub1_clin_net;
sub2_clin_net_thr_d = sub2_clin_net;
str_sub1 = [];
str_sub2 = [];
for i=1:length(density_thr_1)
    sub1_clin_net_thr_d = sub1_clin_net;
    sub2_clin_net_thr_d = sub2_clin_net;
    
    sub1_clin_net_thr_d = threshold_proportional(sub1_clin_net_thr_d,density_thr_1(i));
    sub2_clin_net_thr_d = threshold_proportional(sub2_clin_net_thr_d,density_thr_1(i));
    %sub1_clin_net_thr_d(sub1_clin_net_thr_d<density_thr(i)) = 0;
    %sub2_clin_net_thr_d(sub2_clin_net_thr_d<density_thr(i)) = 0;
    str_sub1(i) = mean(strengths_und(sub1_clin_net_thr_d));
    str_sub2(i) = mean(strengths_und(sub2_clin_net_thr_d));
end

%% plots
h = figure('position',[200,200,400,300])

plot([1:11],str_sub1,'o-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','#D95319','LineWidth',2,'Color','k')
hold on 
plot([1:11],str_sub2,'s-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','#0072BD','LineWidth',2,'Color','k')
xlim([0,12])
ylim([1,10])
set(gca,'xtick',[1:11],'xticklabel',density_thr_1,'FontWeight','bold','FontSize',10)
xlabel('Density Thresholds','FontSize',14)
ylabel('Network Strength','FontSize',14)
legend('Subtype 1','Subtype 2')

[h,p,c,s] = ttest(str_sub1,str_sub2)
s.tstat
p