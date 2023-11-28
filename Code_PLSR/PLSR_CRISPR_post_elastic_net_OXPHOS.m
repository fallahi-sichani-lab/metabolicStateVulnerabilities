clear all;
clc;
close all;

%% Import data
test_data = readtable('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/metabolic_pathways_test_dataset.csv');
train_data = readtable('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/metabolic_pathways_training_dataset.csv');


%% Organize data and create arrays
% converting pathways into categorical data
train_data.Pathway = categorical(train_data.Pathway);

idx_train_OXPHOS = train_data.Pathway == 'Oxidative phosphorylation';


% creating tables based on metabolic state idexesex
OXPHOS_train = train_data(idx_train_OXPHOS,:);


OXPHOS_train_scores = OXPHOS_train(:,1);
OXPHOS_train_Chronos = OXPHOS_train(:,4:end);

% converting pathways into categorical data
test_data.Pathway = categorical(test_data.Pathway);

idx_test_OXPHOS = test_data.Pathway == 'Oxidative phosphorylation';


% creating tables based on metabolic state idex
OXPHOS_test = test_data(idx_test_OXPHOS,:);
OXPHOS_test_scores = OXPHOS_test(:,1);
OXPHOS_test_Chronos = OXPHOS_test(:,4:end);
%% Extracting variables for OXPHOS PLSR
X = table2array(OXPHOS_train_Chronos);
Y = table2array(OXPHOS_train_scores);

%% Model OXPHOS, leave one out CV
clear XLoading YLoading XScores YScores BETA PCTVAR MSE stats; 
ncomp = min(size(X,1)-1,size(X,2));
TSS = sum((Y-mean(Y)).^2);
[XLoading,YLoading,XScores,YScores,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp,'cv',length(Y));
%%
%calculating cummulative variance and varaince across x andy 
var_x0 = PCTVAR(1,:);
cumsum_x0 = cumsum(var_x0)*100;
var_y0 = PCTVAR(2,:);
cumsum_y0 = cumsum(var_y0)*100;

% Prediction accuracy (10-fold cross validation)
Qsquare = [0 1-length(Y)*MSE(2,2:end)/TSS];
% Performance
Rsquare = [0 cumsum(PCTVAR(2,:))];

close all;
f = figure(1); % R2 and Q2 evaluation
hold on
f.Position = [200 200 700 500];


plot(0:ncomp,100*Rsquare,'-b','linewidth',2,'marker','o','markersize',10');
plot(0:ncomp,100*Qsquare,'-r','linewidth',2,'marker','o','markersize',10');
set(gca,'YLim',[0 100],'XLim', [0 20], 'Box','off','XTick',1:2:length(Y));
xlabel('PLS Component');
ylabel('% Variance Explained/Predicted');
legend({'R^2' 'Q^2'});
set(gca,'fontsize',20);

filename = sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/R2_Q2_plot.pdf');
saveas(gcf,filename)

%%
close all;
f = figure(2); 
hold on
f.Position = [200 200 1200 1000];

% PLS1 vs PLS2-PLS5
for i = 2:5
    subplot(2,2,i-1)
    scatter(XScores(:,1), XScores(:,i),200, table2array(OXPHOS_train_scores), 'filled', 'MarkerEdgeColor', 'black');
    hold on;
    
    colorbar
    caxis([-1 1])
    colormap redblue
    xlim([-0.2 0.2])
    ylim([-0.2 0.2])
    xlabel(strcat(strcat("PLS1 (", string(cumsum_y0(1)))," %)"))
    ylabel(strcat(strcat("PLS", string(i)," (", string(cumsum_y0(i)-cumsum_y0(i-1)))," %)"))
    title('PLSR OXPHOS model scores')
    set(gca,'FontSize',14)
end

filename = sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/PLS1_PLS5_scores_plot.pdf');
saveas(gcf,filename)

%% Optimized Model OXPHOS
clear XLoading YLoading XScores YScores BETA PCTVAR MSE stats; 
ncomp = 4;
[XLoading,YLoading,XScores,YScores,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp,'cv',length(Y));
%% Calculate VIP (Variable Importance in Projection) Scores (melanoerentiated model)

[n_observations n_variables] = size(X);

ncomp = 4;

sum1 = zeros(1,n_variables);
sum2 = 0;
clear SS Wnorm2;
for i = 1:ncomp
    SS(i) = (YLoading(i)^2)*(XScores(:,i)'*XScores(:,i));
end
for i = 1:ncomp
    sum2 = sum2 + SS(i);
    Wnorm2(i) = stats.W(:,i)'*stats.W(:,i);
end

clear VIP;
for counter = 1:n_variables
    for k = 1:ncomp
        sum1(counter) = sum1(counter) + SS(k)*stats.W(counter,k)^2/Wnorm2(k);
    end
    VIP(counter) = (n_variables*sum1(counter)/sum2)^0.5;
end

% Plot VIP scores
figure(3);
bar(VIP);
yline(1)
set(gca,'XTick',1:size(X,2),'XTickLabel',OXPHOS_train_Chronos.Properties.VariableNames(:),'Fontsize',12);
title('Variable Importance in Projection','Fontsize',12);
xlabel('Genes','Fontsize',12);
ylabel('VIP score','Fontsize',12);
xtickangle(45);
%% Calculating pearson correlation between metabolic state scores and Chronos scores to calculate signed VIPs

for i = 1:width(OXPHOS_train_Chronos)
    [R,P] = corrcoef(OXPHOS_train.mean_zscore,table2array(OXPHOS_train_Chronos(:,i)));
    R_OXPHOS(i,1) = R(1,2);
    P_OXPHOS(i,1) = P(1,2);
end

%multipling VIPs by the sign of the correlation
Pearson_corr = R_OXPHOS;

Pearson_corr(Pearson_corr<0) = -1;
Pearson_corr(Pearson_corr>0) = 1;

for k = 1:width(OXPHOS_train_Chronos)
    VIP_signed(k,1) = Pearson_corr(k,:)*VIP(:,k);
end

%table of pearson correlation for each gene
Pearson_corr = R_OXPHOS;
OXPHOS_genes = OXPHOS_train.Properties.VariableNames(4:end)';
Pearson_corr_OXPHOS = table(Pearson_corr,P_OXPHOS,OXPHOS_genes, VIP_signed);

%%
%plotting signed VIPs for OXPHOS model
color = [0.26 0.43 0.85; 
    0.90 0.90 0.90; 
    0.79 0.05 0.18];

Pearson_corr_OXPHOS = sortrows(Pearson_corr_OXPHOS, "VIP_signed");

variable_names = Pearson_corr_OXPHOS.OXPHOS_genes(:);

f=figure(4);
hold on;

for k = 1:length(variable_names)
    if Pearson_corr_OXPHOS.VIP_signed(k) < -3
        i=3;
    elseif Pearson_corr_OXPHOS.VIP_signed(k) > 3
        i=1;
    else
        i=2;
    end
    
    h=bar(k,Pearson_corr_OXPHOS.VIP_signed(k));
    set(h,'FaceColor',color(i,:));
end
yline([-0.75 0.75], '-k','LineWidth',1)
yline(1, '-b','LineWidth',1)
yline(-1, '-r','LineWidth',1)
set(gca,'XTick',1:length(variable_names),'XTickLabel','','Fontsize',8);
title('Variable Importance in Projection','Fontsize',20);
xlabel('Gene KO','Fontsize',14);
ylabel('VIP score (signed)','Fontsize',14);
%set(gca, 'FontSize', 8)
xtickangle(45);

hold off;

%%
%exporting table with VIP and Pearson's R
writetable(Pearson_corr_OXPHOS, "OXPHOS_PLSR_selected_variables.txt", 'Delimiter','\t');
%%
%plotting signed VIPs for OXPHOS model
color = [0.26 0.43 0.85; 
    0.90 0.90 0.90; 
    0.79 0.05 0.18];

filtered_VIP = Pearson_corr_OXPHOS(abs(Pearson_corr_OXPHOS.VIP_signed) >= 1,:);
filtered_VIP = sortrows(filtered_VIP, "VIP_signed");

variable_names = filtered_VIP.OXPHOS_genes(:);

f = figure(5);
hold on;
f.Position = [200 200 1000 500];

for k = 1:length(variable_names)
    if filtered_VIP.VIP_signed(k) < -1
        i=1;
    elseif filtered_VIP.VIP_signed(k) > 1
        i=3;
    else
        i=3;
    end
    
    h=bar(k,filtered_VIP.VIP_signed(k));
    set(h,'FaceColor',color(i,:));
end

set(gca,'XTick',1:length(variable_names),'XTickLabel',variable_names(:),'Fontsize',8);
title('Variable Importance in Projection','Fontsize',20);
xlabel('Gene KO','Fontsize',14);
ylabel('VIP score (signed)','Fontsize',14);
xtickangle(45);

filename = sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/VIP_cutoff_1_plot.pdf');
saveas(gcf,filename)

%%
%cross validation of model
[n,p] = size(OXPHOS_train);
yfitPLS = [ones(n,1) table2array(OXPHOS_train_Chronos)]*BETA;

close all;

f = figure(6);
hold on;
f.Position = [200 200 600 500];

mdl = fitlm(OXPHOS_train_scores.mean_zscore(:),yfitPLS);

[R,P] = corrcoef(OXPHOS_train_scores.mean_zscore(:),yfitPLS);
plot([-4 2],[-4 2],'linestyle','-','linewidth',1,'color','k');
hold on;
plot(OXPHOS_train_scores.mean_zscore(:),yfitPLS, 'linestyle','none','marker','o','MarkerSize',16, 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','black') 
xlim([-4 1.2])
ylim([-4 1.2])
set(gca,'fontsize',14);


formatSpec = "PLSR: LOO CV R = %0.2f P = %0.4f";
title(sprintf(formatSpec, R(1,2), P(1,2)), 'Fontsize', 12);
xlabel('Actual OXPHOS score');
ylabel('Predicted OXPHOS score');


filename = sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/LOO_CV_OXPHOS_pancancer_model.pdf');
saveas(gcf,filename)
%%
%independent validation of model
[n,p] = size(OXPHOS_test);
yfitPLS = [ones(n,1) table2array(OXPHOS_test_Chronos)]*BETA;

f = figure(7);
hold on;
f.Position = [200 200 600 500];

mdl = fitlm(OXPHOS_test_scores.mean_zscore(:),yfitPLS);

[R,P] = corrcoef(OXPHOS_test_scores.mean_zscore(:),yfitPLS);
plot([-4 2],[-4 2],'linestyle','-','linewidth',1,'color','k');
hold on;
plot(OXPHOS_test_scores.mean_zscore(:),yfitPLS, 'linestyle','none','marker','o','MarkerSize',20, 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','black') 
set(gca,'fontsize',14);
xlim([-2 1.1])
ylim([-2 1.1])

formatSpec = "PLSR independent validation R = %0.2f P = %0.4f";
title(sprintf(formatSpec, R(1,2), P(1,2)), 'Fontsize', 12);
xlabel('Actual OXPHOS score');
ylabel('Predicted OXPHOS score');

filename = sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/independent_validation_OXPHOS_pancancer_model.pdf');
saveas(gcf,filename)
%% 