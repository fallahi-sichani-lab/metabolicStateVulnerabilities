


clear;
close all;
clc;

drug_sensitivity_AUC_data = readtable('all_cell_lines_Drug_sensitivity_AUC_(CTD^2).xlsx','ReadVariableNames',true);
selected_mutations_data = readtable('all_cell_lines_selected_mutations.xlsx','ReadVariableNames',true);
selected_ssGSE_data = readtable('all_cell_lines_ssGSEA.xlsx','ReadVariableNames',true);


drug_list = (drug_sensitivity_AUC_data(1,9:end).Properties.VariableNames)';
cell_line_list = drug_sensitivity_AUC_data{1:end,2};
cell_line_lineage_list = drug_sensitivity_AUC_data{1:end,4};
AUC = drug_sensitivity_AUC_data{1:end,9:end};
mutation_list = (selected_mutations_data(1,9:end).Properties.VariableNames)';
mutations = selected_mutations_data{1:end,9:end};

ssGSE_list = (selected_ssGSE_data(1,9:end).Properties.VariableNames)';
ssGSE = selected_ssGSE_data{1:end,9:end};

% Identify cell lines included in CTRP dataset and use those in the rest of the analysis
CTRP_cell_line_IDs = (sum((~isnan(AUC))') > 0)';
CTRP_cell_line_list = cell_line_list(CTRP_cell_line_IDs);
CTRP_cell_line_lineage_list = cell_line_lineage_list(CTRP_cell_line_IDs);
CTRP_AUC = AUC(CTRP_cell_line_IDs,:);
CTRP_mutations = mutations(CTRP_cell_line_IDs,:);
CTRP_ssGSE = ssGSE(CTRP_cell_line_IDs,:);


oxphos_high_threshold = prctile(CTRP_ssGSE(:,1),67);
oxphos_low_threshold = prctile(CTRP_ssGSE(:,1),33);


oxphos_high_cell_line_IDs = CTRP_ssGSE(:,1) > oxphos_high_threshold;
oxphos_high_cell_line_list = CTRP_cell_line_list(oxphos_high_cell_line_IDs);
oxphos_low_cell_line_IDs = CTRP_ssGSE(:,1) < oxphos_low_threshold;
oxphos_low_cell_line_list = CTRP_cell_line_list(oxphos_low_cell_line_IDs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identifying drugs with enhanced effect on Oxphos high lines in the presence of PTEN mutations

PTEN_damaging_mutated_cell_line_IDs = CTRP_mutations(:,6) > 0;
PTEN_damaging_mutated_cell_lines = CTRP_cell_line_list(PTEN_damaging_mutated_cell_line_IDs);

PTEN_WT_cell_line_IDs = (CTRP_mutations(:,6) + CTRP_mutations(:,11) == 0);
PTEN_WT_cell_lines = CTRP_cell_line_list(PTEN_WT_cell_line_IDs);

oxphos_high_PTEN_mutated_cell_line_IDs = oxphos_high_cell_line_IDs .* PTEN_damaging_mutated_cell_line_IDs;
oxphos_high_PTEN_mutated_cell_lines = CTRP_cell_line_list(oxphos_high_PTEN_mutated_cell_line_IDs > 0);

oxphos_high_PTEN_WT_cell_line_IDs = oxphos_high_cell_line_IDs .* PTEN_WT_cell_line_IDs;
oxphos_high_PTEN_WT_cell_lines = CTRP_cell_line_list(oxphos_high_PTEN_WT_cell_line_IDs > 0);


clear AUC_oxphos_low_cell_lines AUC_oxphos_high_cell_lines AUC_PTEN_damaging_mutated_cell_lines AUC_oxphos_high_PTEN_mutated_cell_lines AUC_oxphos_high_PTEN_WT_cell_lines;
for drug_id = 1:length(drug_list)
    AUC_oxphos_low_cell_lines(drug_id,1:length(oxphos_low_cell_line_list)) = CTRP_AUC(oxphos_low_cell_line_IDs,drug_id)';
    AUC_oxphos_high_cell_lines(drug_id,1:length(oxphos_high_cell_line_list)) = CTRP_AUC(oxphos_high_cell_line_IDs,drug_id)';
    AUC_oxphos_high_PTEN_mutated_cell_lines(drug_id,1:length(oxphos_high_PTEN_mutated_cell_lines)) = CTRP_AUC(oxphos_high_PTEN_mutated_cell_line_IDs > 0,drug_id)';
    AUC_oxphos_high_PTEN_WT_cell_lines(drug_id,1:length(oxphos_high_PTEN_WT_cell_lines)) = CTRP_AUC(oxphos_high_PTEN_WT_cell_line_IDs > 0,drug_id)';
end

figure(1);
plot([0 0],[-3 3],'--k','LineWidth',1);
hold on;
plot([-2 2],[0 0],'--k','LineWidth',1);

clear differential_AUC_oxphos_low_vs_high differential_AUC_oxphos_high_PTEN_WT_vs_MUT;
clear differential_AUC_oxphos_low_vs_high_pvalue differential_AUC_oxphos_high_PTEN_WT_vs_MUT_pvalue;
oxphos_high_PTEN_mutation_synthetic_lethal_drug_IDs = [];

for drug_id = 1:length(drug_list)
    differential_AUC_oxphos_low_vs_high(drug_id) = nanmedian(AUC_oxphos_low_cell_lines(drug_id,:)) - nanmedian(AUC_oxphos_high_cell_lines(drug_id,:));
    [p h] = ranksum(AUC_oxphos_low_cell_lines(drug_id,:),AUC_oxphos_high_cell_lines(drug_id,:),'tail','right');
    differential_AUC_oxphos_low_vs_high_pvalue(drug_id) = p;

    differential_AUC_oxphos_high_PTEN_WT_vs_MUT(drug_id) = nanmedian(AUC_oxphos_high_PTEN_WT_cell_lines(drug_id,:)) - nanmedian(AUC_oxphos_high_PTEN_mutated_cell_lines(drug_id,:));
    [p h] = ranksum(AUC_oxphos_high_PTEN_WT_cell_lines(drug_id,:),AUC_oxphos_high_PTEN_mutated_cell_lines(drug_id,:),'tail','right');
    differential_AUC_oxphos_high_PTEN_WT_vs_MUT_pvalue(drug_id) = p;
end

for drug_id = 1:length(drug_list)
    if ((differential_AUC_oxphos_low_vs_high(drug_id) > 0.5) && (differential_AUC_oxphos_high_PTEN_WT_vs_MUT(drug_id) > 0.5) &&...
            (differential_AUC_oxphos_low_vs_high_pvalue(drug_id) < 0.05))
        oxphos_high_PTEN_mutation_synthetic_lethal_drug_IDs = [oxphos_high_PTEN_mutation_synthetic_lethal_drug_IDs; drug_id];
        h = plot(differential_AUC_oxphos_low_vs_high(drug_id),differential_AUC_oxphos_high_PTEN_WT_vs_MUT(drug_id),...
            'marker','o','markersize',12,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0]);
%         gname(drug_list(drug_id));
    else
        h = plot(differential_AUC_oxphos_low_vs_high(drug_id),differential_AUC_oxphos_high_PTEN_WT_vs_MUT(drug_id),...
            'marker','o','markersize',7,'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0 0 0]);
    end
end
set(gca,'Box','off','FontSize',15,'XLim',[-2 2],'YLim',[-3 3]);
title ('Differential median AUC')
xlabel('OxhphosLow - OxhphosHigh');
ylabel('OxhphosHigh: PTEN(WT) - PTEN(MUT)');

drug_list(oxphos_high_PTEN_mutation_synthetic_lethal_drug_IDs)


