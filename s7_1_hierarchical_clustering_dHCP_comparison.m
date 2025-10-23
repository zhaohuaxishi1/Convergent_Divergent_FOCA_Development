%% ================================================================
% Script: hierarchical_clustering_dHCP_comparison.m
% Purpose:
%   Compare neonatal FOCA hierarchical structure with adult-derived
%   clustering organization (MSC reference).
%
% Description:
%   - Compute neonatal FOCA matrix from individualized network topographies.
%   - Apply adult hierarchical clustering solution (Z, outperm) for alignment.
%   - Visualize reordered FOCA matrix and summarize macro-scale modularity.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
%
% ================================================================

clc; clear;
addpath(genpath('./dependence/customcolormap'));

%% --- Step 1. Load neonatal individualized topography ---
disp('--- Loading neonatal group-level topography ---');
load('./data/dHCP/uCi_Group_Topography_dHCP');
data = uCi_Group_Topography_dHCP;
nNet = size(data, 2);

%% --- Step 2. Compute FOCA matrix for neonates ---
disp('--- Computing FOCA (network × network correlation matrix) ---');
[corrMatrix_dHCP, pval_dHCP] = corr(data);
distMatrix_dHCP = 1 - corrMatrix_dHCP;
save('./results/FOCA/FOCA_GroupLevel_dHCP.mat', 'corrMatrix_dHCP', 'pval_dHCP', 'distMatrix_dHCP');

%% --- Step 3. Load adult (MSC) hierarchical clustering structure ---
disp('--- Loading adult (MSC) hierarchical clustering structure ---');

% Fallback check
if ~exist('Z_MSC','var')
    load('./results/FOCA/MSC_FOCA_Hierarchical_Structure.mat', 'Z', 'outperm');
    Z_MSC = Z;
    outperm_MSC = outperm;
end

%% --- Step 4. Apply adult-derived cluster order to neonatal FOCA ---
disp('--- Applying adult-derived ordering to neonatal FOCA ---');
corrMatrix_dHCP_reordered = corrMatrix_dHCP(outperm_MSC, outperm_MSC);

%% --- Step 5. Visualize reordered neonatal FOCA matrix ---
disp('--- Visualizing reordered neonatal FOCA matrix ---');
figure('Name','Neonatal FOCA (Adult Cluster Order)',...
    'NumberTitle','off','Units','inches','Position',[0 0 9 8]);
imagesc(corrMatrix_dHCP_reordered);
colormap(customcolormap_preset('red-white-blue'));
caxis([-0.6 1]);
cb = colorbar; cb.Ticks = -0.6:0.4:1;

network_labels = {'DMN Par','DMN Ant','DMN Dor','DMN Ret',...
                  'Vis Lat','Vis Str','Vis V5','Vis V1',...
                  'FP','DAN','PreMot','Lang','Sal','CO/AMN',...
                  'PMN','Hand','Face','Foot','Aud','SCAN'};
newLabels = network_labels(outperm_MSC);
set(gca,'XTick',1:nNet,'YTick',1:nNet,...
    'XTickLabel',newLabels,'YTickLabel',newLabels,...
    'FontName','Arial','FontSize',14,'TickLength',[0 0]);
xtickangle(90); axis square;
title('Neonatal FOCA Matrix (Adult Cluster Order)');
print(gcf, './results/FOCA/FOCA_ReorderedMatrix_dHCP.png', '-r300', '-dpng');
close(gcf);

%% --- Step 6. Summarize large-scale modular organization ---
disp('--- Summarizing macro-scale modular organization ---');
m = corrMatrix_dHCP_reordered;

% Use adult-defined macro-cluster mapping (for developmental comparison)
corrMatrixLabel_MSC = [1,1,1,1,2,2,3,3,3,3,3,4,4,4,4,5,5,5,6,6];
uCi_macro = unique(nonzeros(corrMatrixLabel_MSC));
uCi_FC_dHCP = zeros(length(uCi_macro));

for i = 1:length(uCi_macro)
    for j = 1:length(uCi_macro)
        uCi_FC_dHCP(i,j) = nanmean(nanmean(m(corrMatrixLabel_MSC==uCi_macro(i), corrMatrixLabel_MSC==uCi_macro(j))));
    end
end

macroLabels = {'DMN','PMN-Sal','Vis-Lang','Mot','AMN','FP-DAN'};
save('./results/FOCA/FOCA_Modular_Summary_dHCP.mat', 'uCi_FC_dHCP', 'macroLabels');

figure('Name','Neonatal FOCA Macro Modularity',...
    'NumberTitle','off','Units','inches','Position',[0 0 6 5]);
imagesc(uCi_FC_dHCP);
colormap(customcolormap_preset('red-white-blue'));
caxis([-0.3 0.5]);
cb = colorbar; cb.Ticks = -0.3:0.2:0.5;
set(gca,'XTick',1:length(macroLabels),'XTickLabel',macroLabels,...
    'YTick',1:length(macroLabels),'YTickLabel',macroLabels,...
    'FontName','Arial','FontSize',14,'TickLength',[0 0]);
xtickangle(90); axis square;
title('Neonatal Modular Organization (Adult-defined Systems)');
print(gcf, './results/FOCA/FOCA_Modular_Summary_dHCP.png', '-r300', '-dpng');

disp('--- Neonatal FOCA hierarchical clustering comparison complete ---');
