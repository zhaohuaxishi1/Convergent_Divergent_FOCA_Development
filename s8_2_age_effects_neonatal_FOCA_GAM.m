%% ================================================================
% Script: age_effects_neonatal_FOCA_GAM.m
% Purpose:
%   Model age-related trajectories of internetwork covariance (FOCA matrix)
%   during the early postnatal period using mass univariate GAMs.
%
% Description:
%   - Input: precomputed F-statistic matrix (Gam_F_Vector_all, 20×20)
%            and p-value matrix (Gam_P_Vector_all, 20×20)
%   - Apply Benjamini–Hochberg FDR correction (q < 0.05)
%   - Generate significance mask and symmetrize matrices
%   - Visualize F-statistics reordered by adult (MSC) hierarchical clustering
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;
addpath(genpath('./dependence/customcolormap'));

%% --- Step 1. Load precomputed GAM results (F-statistics) ---
disp('--- Loading precomputed GAM results ---');
load('./data/dHCP/GAM_AgeEffects/Gam_F_Vector_all.mat', 'Gam_F_Vector_all');

if ~exist('Gam_F_Vector_all','var')
    error('Variable Gam_F_Vector_all not found.');
end

nNet = size(Gam_F_Vector_all,1);
%% --- Step 3. Reorder using adult (MSC) hierarchical structure ---
disp('--- Reordering using adult (MSC) FOCA hierarchy ---');
% Fallback check
if ~exist('Z_MSC','var')
    load('./results/FOCA/MSC_FOCA_Hierarchical_Structure.mat', 'Z', 'outperm');
    Z_MSC = Z;
    outperm_MSC = outperm;
end

Gam_F_Reordered = Gam_F_Vector_all(outperm_MSC, outperm_MSC);

%% --- Step 4. Visualization ---
disp('--- Visualizing neonatal FOCA age-related effects (GAM F-statistics) ---');
network_labels = {'DMN Par','DMN Ant','DMN Dor','DMN Ret', ...
                  'Vis Lat','Vis Str','Vis V5','Vis V1', ...
                  'FP','DAN','PreMot','Lang','Sal','CO/AMN', ...
                  'PMN','Hand','Face','Foot','Aud','SCAN'};
newLabels = network_labels(outperm_MSC);

figure('Name','Neonatal FOCA Age Effects (GAM F-statistics)', ...
    'NumberTitle','off','Units','inches','Position',[0 0 9 8]);
imagesc(Gam_F_Reordered);
myCmap = colormap(customcolormap_preset('orange-white-purple'));
myCmap(1:130,:) = [];
colormap(myCmap);

caxis([0 10]);
cb = colorbar;
cb.Ticks = 0:2:10;
cb.Label.String = 'F-statistics (GAM)';
set(gca,'XTick',1:nNet,'YTick',1:nNet,...
    'XTickLabel',newLabels,'YTickLabel',newLabels,...
    'FontName','Arial','FontSize',18,'TickLength',[0 0]);
xtickangle(90);
axis square;
title('Age-related Effects on FOCA (Neonates, GAM)');

%% --- Step 5. Save results ---
outputDir = './results/GAM/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

print(gcf, [outputDir 'dHCP_FOCA_AgeEffects_GAM_F_Reordered.png'], '-r300', '-dpng');

close(gcf);
disp('--- GAM age-effect visualization complete ---');
