%% cell type analysis
clear
clc

load('100DS82scaledRobustSigmoidNSGRNAseqQC1Lcortex_ROI_NOdistCorrEuclidean(1).mat');
load('cell_genes_5.mat')

[~,all_genes] = xlsread('all_genes.xlsx');

load('pls_genes.mat')
ex_n(762)=[];
in_n(201)=[];
genes = probeInformation.GeneSymbol;

disp('finished......')
%% calculation cell type enrichments

[cell_gene_num_sub1,p_permu_sub1,p_pdr_sub1,cell_type_genes_sub1] = cell_type_enrichment(genes,sub1_pls_genes,astro,endo,ex_n,in_n,micro,olig,opc);
[cell_gene_num_sub2,p_permu_sub2,p_pdr_sub2,cell_type_genes_sub2] = cell_type_enrichment(genes,sub2_pls_genes,astro,endo,ex_n,in_n,micro,olig,opc);

disp('finished......')
%% write overlapped genes cell types 

cell_types = {'Astro','Endo','Neuro.ex','Neuro.in','Micro','Oligo','OPC'};

fid1 = fopen(['subtype1_cell_type_overlapped_genes.csv'],'w');
fprintf(fid1,'%s,%s\n','Cell Types','Gene names');
for i=1:7
    genes_name = ',';
    for j=1:length(cell_type_genes_sub1{i})
       genes_name = strcat(genes_name,',',cell_type_genes_sub1{i}{j});
    end
    fprintf(fid1,'%s,%s\n',cell_types{i},genes_name);
end
fclose(fid1);

fid1 = fopen(['subtype2_cell_type_overlapped_genes.csv'],'w');
fprintf(fid1,'%s,%s\n','Cell Types','Gene names');
for i=1:7
    genes_name = ',';
    for j=1:length(cell_type_genes_sub2{i})
       genes_name = strcat(genes_name,',',cell_type_genes_sub2{i}{j});
    end
    fprintf(fid1,'%s,%s\n',cell_types{i},genes_name);
end
fclose(fid1);

disp('finished......')
%% overlap ratios  

[permut_over]=permut_rate(genes,sub2_pls_genes,micro);

permut_over = permut_over';
subtype2_in_permut5000 = permut_over'/cell_gene_num_sub2(2,4);

cell_gene_num_sub1(3,:) = cell_gene_num_sub1(1,:)./cell_gene_num_sub1(2,:);
cell_gene_num_sub2(3,:) = cell_gene_num_sub2(1,:)./cell_gene_num_sub2(2,:);

disp('finished......')
%% chi square test between subgroup 1 and subgroup 2
%%%% Table S3

chi_f=[];
chi_p = [];

for i=1:7
    a = [length(sub1_pls_genes)-cell_gene_num_sub1(1,i),cell_gene_num_sub1(1,i);length(sub2_pls_genes)-cell_gene_num_sub2(1,i),cell_gene_num_sub2(1,i)];
    %a = [cell_gene_num_sub1(1,i),length(sub1_pls_genes)-cell_gene_num_sub1(1,i);cell_gene_num_sub1(2,i)-cell_gene_num_sub1(1,i),10027-length(sub1_pls_genes)-cell_gene_num_sub1(2,i)+cell_gene_num_sub1(1,i)];
    [chi_p(i), chi_f(i)]= chi2test(a);
end

chi_p = mafdr(chi_p,'BHFDR',true);
chi_p = chi_p';
disp('finished......')
%%
% chi_f1=[];
% chi_p1 = [];
% 
% for i=1:7
%     %a = [length(sub1_pls_genes),cell_gene_num_sub1(1,i);length(sub2_pls_genes),cell_gene_num_sub2(1,i)];
%     a = [cell_gene_num_sub1(1,i),length(sub1_pls_genes)-cell_gene_num_sub1(1,i);cell_gene_num_sub1(2,i)-cell_gene_num_sub1(1,i),10027-length(sub1_pls_genes)-cell_gene_num_sub1(2,i)+cell_gene_num_sub1(1,i)];
%     %[chi_p1(i), Q1]= chi2test(a);
%     [h,chi_p1(i),stats] = fishertest(a);
% end
% 
% chi_f2=[];
% chi_p2 = [];
% 
% for i=1:7
%     %a = [length(sub1_pls_genes),cell_gene_num_sub1(1,i);length(sub2_pls_genes),cell_gene_num_sub2(1,i)];
%     a = [cell_gene_num_sub2(1,i),length(sub2_pls_genes)-cell_gene_num_sub2(1,i);cell_gene_num_sub2(2,i)-cell_gene_num_sub2(1,i),10027-length(sub2_pls_genes)-cell_gene_num_sub2(2,i)+cell_gene_num_sub2(1,i)];
%     %[chi_p2(i), Q2]= chi2test(a);
%     [h,chi_p2(i),stats] = fishertest(a);
% end

%%
cell_gene_num_sub1(4,:) = cell_gene_num_sub1(1,:)./length(sub1_pls_genes);
cell_gene_num_sub2(4,:) = cell_gene_num_sub2(1,:)./length(sub2_pls_genes);