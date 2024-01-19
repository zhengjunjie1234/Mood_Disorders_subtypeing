function [cell_gene_num,p_permu,pfdr,cell_type_genes] = cell_type_enrichment(genes,gene_list,astro,endo,ex_n,in_n,micro,olig,opc)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%  genes = backgroud genes 
%  genelist = interest genes
%
%%% permutation test
astro_u = intersect(genes,astro);
%astro_u = astro;
%[p_permu_astro,astro_in] = permutation_over(genes,astro_u,gene_list);
[p_permu_astro,astro_in] = permutation_over_new(genes,astro_u,gene_list);

endo_u = intersect(genes,endo);
%endo_u = endo;
%[p_permu_endo,endo_in] = permutation_over(genes,endo_u,gene_list);
[p_permu_endo,endo_in] = permutation_over_new(genes,endo_u,gene_list);

ex_n_u = intersect(genes,ex_n);
%ex_n_u = ex_n;
%[p_permu_ex_n,ex_in] = permutation_over(genes,ex_n_u,gene_list);
[p_permu_ex_n,ex_in] = permutation_over_new(genes,ex_n_u,gene_list);

in_n_u = intersect(genes,in_n);
%in_n_u = in_n;
%[p_permu_in_n,in_in] = permutation_over(genes,in_n_u,gene_list);
[p_permu_in_n,in_in] = permutation_over_new(genes,in_n_u,gene_list);

micro_u =intersect(genes,micro);
%micro_u = micro;
%[p_permu_micro,micro_in] = permutation_over(genes,micro_u,gene_list);
[p_permu_micro,micro_in] = permutation_over_new(genes,micro_u,gene_list);

olig_u = intersect(genes,olig);
%olig_u = olig;
%[p_permu_olig,olig_in] = permutation_over(genes,olig_u,gene_list);
[p_permu_olig,olig_in] = permutation_over_new(genes,olig_u,gene_list);

opc_u = intersect(genes,opc);
%opc_u = opc;
%[p_permu_opc_u,opc_in] = permutation_over(genes,opc_u,gene_list);
[p_permu_opc_u,opc_in] = permutation_over_new(genes,opc_u,gene_list);

cell_gene_num(1,:) = [length(astro_in),length(endo_in),length(ex_in),length(in_in),length(micro_in),length(olig_in),length(opc_in)];
cell_gene_num(2,:) = [length(astro_u),length(endo_u),length(ex_n_u),length(in_n_u),length(micro_u),length(olig_u),length(opc_u)];


pfdr = mafdr([p_permu_astro,p_permu_endo,p_permu_ex_n,p_permu_in_n,p_permu_micro,p_permu_olig,p_permu_opc_u],'BHFDR',true);
p_permu = [p_permu_astro,p_permu_endo,p_permu_ex_n,p_permu_in_n,p_permu_micro,p_permu_olig,p_permu_opc_u];

cell_type_genes{1} = astro_in;
cell_type_genes{2} = endo_in;
cell_type_genes{3} = ex_in;
cell_type_genes{4} = in_in;
cell_type_genes{5} = micro_in;
cell_type_genes{6} = olig_in;
cell_type_genes{7} = opc_in;

end

