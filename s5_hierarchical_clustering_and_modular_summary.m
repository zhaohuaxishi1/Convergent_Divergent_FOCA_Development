%% ================================================================
% Script: hierarchical_clustering_and_modular_summary.m
% Purpose:
%   Quantify hierarchical and modular organization of internetwork 
%   topographic covariance (FOCA matrix) in the MSC dataset.
%
% Description:
%   - Compute pairwise correlations (FOCA matrix) among 20 networks.
%   - Perform hierarchical clustering using average linkage (1 - r).
%   - Visualize dendrogram and reordered FOCA matrix.
%   - Summarize intra- and inter-cluster correlations across large-scale systems.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
%
% Note:
%   Clustering was performed using Matlab R2020b, average linkage,
%   and distance metric (1 - Pearson correlation).
% ================================================================

clc; clear;
addpath(genpath('./dependence/customcolormap'));

%% --- Step 1. Load and average individualized topographies (MSC dataset) ---
dataFile = './data/';
load([dataFile, 'Individual_Topography_Lynch20_MSC_all.mat']);
nNet = 20;
nVertex = 59412;

disp('--- Averaging individualized functional topographies across subjects ---');
uCi_Group_Topography = zeros(nVertex, nNet);
validCount = 0;

for i = 1:length(uCi_individual_network_all)
    if size(uCi_individual_network_all{i}, 2) == nNet
        uCi_Group_Topography = uCi_Group_Topography + uCi_individual_network_all{i};
        validCount = validCount + 1;
    end
end
uCi_Group_Topography = uCi_Group_Topography ./ validCount;
fprintf('Averaged %d subjects.\n', validCount);

%% --- Step 2. Compute FOCA matrix (network × network correlations) ---
disp('--- Computing FOCA (network-level spatial correlation matrix) ---');
[corrMatrix, pval] = corr(uCi_Group_Topography);
distMatrix = 1 - corrMatrix;
save('./results/FOCA/MSC_FOCA_GroupLevel_Matrix.mat', 'corrMatrix', 'pval', 'distMatrix');

%% --- Step 3. Hierarchical clustering ---
disp('--- Performing hierarchical clustering ---');
Y = squareform(distMatrix, 'tovector');   % convert to condensed distance vector
Z = linkage(Y, 'average');                % average linkage method

% Dendrogram
figure('Name','Hierarchical Clustering Dendrogram','NumberTitle','off','Units','inches','Position',[0 0 8 2]);
[H,~,outperm] = dendrogram(Z, 0, 'Orientation', 'top');
set(H, 'LineWidth', 2);
title('Hierarchical Clustering of FOCA Matrix');
xlabel('Functional Networks'); ylabel('1 - Pearson Correlation');
set(gca, 'FontName','Arial','FontSize',12,'TickDir','out');
save('./results/FOCA/MSC_FOCA_Hierarchical_Structure.mat', 'Z', 'outperm');

%% --- Step 4. Reorder and visualize FOCA matrix ---
disp('--- Visualizing reordered FOCA correlation matrix ---');
corrMatrixReordered = corrMatrix(outperm, outperm);
figure('Name','FOCA Reordered Correlation Matrix','NumberTitle','off','Units','inches','Position',[0 0 9 8]);
imagesc(corrMatrixReordered);
colormap(customcolormap_preset('red-white-blue'));
caxis([-0.6 1]);
cb = colorbar; cb.Ticks = -0.6:0.4:1;
network_labels = {'DMN Par','DMN Ant','DMN Dor','DMN Ret', ...
                  'Vis Lat','Vis Str','Vis V5','Vis V1', ...
                  'FP','DAN','PreMot','Lang','Sal','CO/AMN', ...
                  'PMN','Hand','Face','Foot','Aud','SCAN'};
newLabels = network_labels(outperm);
set(gca,'XTick',1:nNet,'YTick',1:nNet,...
    'XTickLabel',newLabels,'YTickLabel',newLabels,...
    'FontName','Arial','FontSize',12,'TickLength',[0 0]);
xtickangle(90); axis square;
title('Reordered FOCA Correlation Matrix');
print(gcf, './results/FOCA/MSC_FOCA_ReorderedMatrix.png', '-r300', '-dpng');
close(gcf);

%% --- Step 5. Summarize modular organization across large-scale systems ---
disp('--- Summarizing intra- and inter-cluster relationships ---');
m = corrMatrixReordered;

% Define 6 macro-systems based on observed clustering
corrMatrixLabel = [1,1,1,1,2,2,3,3,3,3,3,4,4,4,4,5,5,5,6,6];
uCi = unique(corrMatrixLabel);
uCi_FC = zeros(length(uCi));

for i = 1:length(uCi)
    for j = 1:length(uCi)
        uCi_FC(i,j) = nanmean(nanmean(m(corrMatrixLabel==uCi(i), corrMatrixLabel==uCi(j))));
    end
end

% Plot summary modular FOCA matrix
macroLabels = {'DMN','PMN-Sal','Vis-Lang','Somatomotor','AMN','FP-DAN'};
figure('Name','Macro-level FOCA Modularity','NumberTitle','off','Units','inches','Position',[0 0 6 5]);
imagesc(uCi_FC);
colormap(customcolormap_preset('red-white-blue'));
caxis([-0.3 0.5]); cb = colorbar;
set(gca,'XTick',1:length(macroLabels),'XTickLabel',macroLabels,...
    'YTick',1:length(macroLabels),'YTickLabel',macroLabels,...
    'FontName','Arial','FontSize',14);
xtickangle(90); axis square;
title('Large-scale Modular Organization (FOCA)');
print(gcf, './results/FOCA/MSC_FOCA_Modular_Summary.png', '-r300', '-dpng');

save('./results/FOCA/MSC_FOCA_Modular_Summary.mat', 'uCi_FC', 'macroLabels');
disp('--- Hierarchical clustering and modular summary complete ---');
