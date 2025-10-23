%% ================================================================
% Script: visualize_developmental_difference_FOCA_CohenD.m
% Purpose:
%   Visualize developmental differences in individualized FOCA matrices
%   between neonates (dHCP) and adults (HCP) using precomputed Cohen’s d.
%
% Description:
%   - Load precomputed Cohen’s d matrix (20 × 20).
%   - Reorder using adult hierarchical clustering solution (MSC dataset).
%   - Visualize developmental differences in FOCA using a diverging colormap.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
%
% ================================================================

clc; clear;
addpath(genpath('./dependence/customcolormap'));

%% --- Step 1. Load precomputed results ---
disp('--- Loading precomputed developmental FOCA effect sizes ---');
load('./data/dHCP_HCP_FOCA_CohenD.mat', 'Cohen_d_FNC');  

% Fallback check
if ~exist('Z_MSC','var')
    load('./results/FOCA/MSC_FOCA_Hierarchical_Structure.mat', 'Z', 'outperm');
    Z_MSC = Z;
    outperm_MSC = outperm;
end

if ~exist('Cohen_d_FNC','var')
    error('Cohen_d_FNC variable not found in the loaded file.');
end

nNet = size(Cohen_d_FNC,1);
disp(['Loaded Cohen’s d matrix of size ', num2str(nNet), ' × ', num2str(nNet)]);

%% --- Step 2. Apply adult (MSC) hierarchical reordering ---
disp('--- Reordering matrix according to adult (MSC) hierarchy ---');
CohenD_Reordered = Cohen_d_FNC(outperm_MSC, outperm_MSC);

%% --- Step 3. Visualization ---
disp('--- Visualizing developmental FOCA differences (Cohen’s d) ---');
figure('Name','Developmental FOCA Differences (Cohen’s d)', ...
    'NumberTitle','off','Units','inches','Position',[0 0 9 8]);

imagesc(CohenD_Reordered);
colormap(customcolormap_preset('orange-white-purple'));
caxis([-3 3]);
cb = colorbar;
cb.Ticks = -3:1:3;
cb.Label.String = 'Cohen''s d (dHCP - HCP)';

network_labels = {'DMN Par','DMN Ant','DMN Dor','DMN Ret', ...
                  'Vis Lat','Vis Str','Vis V5','Vis V1', ...
                  'FP','DAN','PreMot','Lang','Sal','CO/AMN', ...
                  'PMN','Hand','Face','Foot','Aud','SCAN'};

newLabels = network_labels(outperm_MSC);
set(gca,'XTick',1:nNet,'YTick',1:nNet,...
    'XTickLabel',newLabels,'YTickLabel',newLabels,...
    'FontName','Arial','FontSize',18,'TickLength',[0 0]);
xtickangle(90); axis square;
title('Developmental Differences in FOCA Matrix');

%% --- Step 4. Save results ---
disp('--- Saving visualization and reordered data ---');
outputDir = './results/FOCA/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

print(gcf, [outputDir 'FOCA_Developmental_Diff_CohenD_Reordered.png'], '-r300', '-dpng');
save([outputDir 'FOCA_Developmental_Diff_CohenD_Reordered.mat'], ...
    'CohenD_Reordered', 'Cohen_d_FNC', 'network_labels', 'outperm_MSC');

close(gcf);
disp('--- Visualization complete ---');
