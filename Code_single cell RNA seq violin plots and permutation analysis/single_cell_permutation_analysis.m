clear all;
clc;
close all;

%% Import data
mean_single_cell_exp = readtable('/Volumes/FallahiLab/Maize-Data/Data/Cara/CellxGene single cell data/mean_mito_OXPHOS_single_cell_exp_dataset.csv');

%% extracting different gene sets/pathways
pathways = unique(mean_single_cell_exp.gene_set);

%%
%permutation test 
for i = 1:length(pathways)

    pathway_idx = strcmp(mean_single_cell_exp.gene_set, pathways(i));
    pathway_subset = mean_single_cell_exp(pathway_idx,:);

    mutated_idx = strcmp(pathway_subset.PTEN_status, 'Mutated');
    mutation_samples = pathway_subset(mutated_idx,:);

    not_mutated_idx = strcmp(pathway_subset.PTEN_status, 'Not Mutated');
    no_mutation_samples = pathway_subset(not_mutated_idx,:);

    [p(i), ~, ~] = permutationTest(mutation_samples.Mean_exp, no_mutation_samples.Mean_exp, 5000, 'sidedness','larger');
end

perm_results = table(p(:),pathways(:),'VariableNames',{'p-value','pathway'});
%%
