%% ================================================================
% Script: FOCA_vs_FC_stability_comparison.m
% Purpose:
%   Compare stability of FOCA (Functional Topographic Covariance Analysis)
%   with conventional FC (Functional Connectivity) under GSR vs. non-GSR.
%
% Description:
%   - Compute mean-normalized RMSE and correlation between GSR/non-GSR maps
%   - Group-level differences tested via two-sample t-test and Cohen's d
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;
addpath(genpath('./dependence'));
% addpath('./utilities');

network_names = { 'DMN Par','DMN Ant','DMN Dor','DMN Ret',...
                  'Vis Lat','Vis Str','Vis V5','Vis V1',...
                  'FP','DAN','PreMot','Lang','Sal','AMN',...
                  'PMN','Hand','Face','Foot','Aud','SCAN'};

mask = triu(true(20,20),1);  % upper triangular mask (exclude diagonal)

%% --- Step 1. Load FOCA data (GSR and non-GSR) ---
disp('--- Loading FOCA data (GSR / non-GSR) ---');

load('./data/HCP/FOCA_Individual_HCP_GSR.mat');     % GSR version
FOCA_Individual_GSR = FOCA_Individual_HCP_GSR;

load('./data/HCP/FOCA_Individual_HCP_NGSR.mat');    % non-GSR version
FOCA_Individual_NGSR = FOCA_Individual_HCP_NGSR;
clear FOCA_Individual_HCP_GSR FOCA_Individual_HCP_NGSR

% Remove empty subjects
valid_idx = find(sum(FOCA_Individual_GSR,2)~=0 & sum(FOCA_Individual_NGSR,2)~=0);
FOCA_Individual_GSR = FOCA_Individual_GSR(valid_idx,:);
FOCA_Individual_NGSR = FOCA_Individual_NGSR(valid_idx,:);

%% --- Step 2. Compute FOCA stability metrics ---
disp('--- Computing FOCA stability metrics ---');
nSub = size(FOCA_Individual_GSR,1);
r_foca = nan(nSub,1);
nrmse_foca = nan(nSub,1);

for s = 1:nSub
    A = reshape(FOCA_Individual_GSR(s,:), 20, 20);
    B = reshape(FOCA_Individual_NGSR(s,:), 20, 20);
    [r_foca(s), ~] = corr(A(mask), B(mask), 'rows','complete');
    diffAB = A(mask) - B(mask);
    RMSE = sqrt(mean(diffAB.^2));
    mean_val = mean(abs([A(mask); B(mask)]));
    nrmse_foca(s) = RMSE / mean_val;
end
fprintf('FOCA: mean r = %.3f, mean nRMSE = %.3f\n', mean(r_foca), mean(nrmse_foca));

%% --- Step 3. Load conventional FC data ---
disp('--- Loading FC data (GSR / non-GSR) ---');

load('./data/HCP/FC_Individual_HCP_GSR.mat');     % GSR version
FC_Individual_GSR = FC_Individual_HCP_GSR;

load('./data/HCP/FC_Individual_HCP_NGSR.mat');    % non-GSR version
FC_Individual_NGSR = FC_Individual_HCP_NGSR;
clear FC_Individual_HCP_GSR FC_Individual_HCP_NGSR

% Remove empty subjects
valid_idx = find(sum(FC_Individual_GSR,2)~=0 & sum(FC_Individual_NGSR,2)~=0);
FC_Individual_GSR = FC_Individual_GSR(valid_idx,:);
FC_Individual_NGSR = FC_Individual_NGSR(valid_idx,:);

%% --- Step 4. Compute FC stability metrics ---
disp('--- Computing FC stability metrics ---');
nSub = size(FC_Individual_GSR,1);
r_fc = nan(nSub,1);
nrmse_fc = nan(nSub,1);

for s = 1:nSub
    A = reshape(FC_Individual_GSR(s,:), 20, 20);
    B = reshape(FC_Individual_NGSR(s,:), 20, 20);
    [r_fc(s), ~] = corr(A(mask), B(mask), 'rows','complete');
    diffAB = A(mask) - B(mask);
    RMSE = sqrt(mean(diffAB.^2));
    mean_val = mean(abs([A(mask); B(mask)]));
    nrmse_fc(s) = RMSE / mean_val;
end
fprintf('FC: mean r = %.3f, mean nRMSE = %.3f\n', mean(r_fc), mean(nrmse_fc));

%% --- Step 5. Group-level comparison (t-test + effect size) ---
disp('--- Group-level comparison: FOCA vs FC stability ---');
[~, p_ttest, ~, stats] = ttest2(nrmse_foca, nrmse_fc);
T_value = stats.tstat;
Cohen_d_value = computeCohen_d(nrmse_foca, nrmse_fc);

fprintf('Two-sample t-test: T = %.3f, p = %.3e, Cohen''s d = %.3f\n', ...
        T_value, p_ttest, Cohen_d_value);

save('./results/StabilityComparison/FOCA_FC_Stability_Stats.mat', ...
     'nrmse_foca', 'nrmse_fc', 'r_foca', 'r_fc', 'T_value', 'p_ttest', 'Cohen_d_value');

%% --- Step 6. Visualization ---
disp('--- Plotting FOCA vs FC normalized RMSE ---');
C = [0 0 0];

figure('Units','centimeters','Position',[0 0 6 7]);
hold on;

pos_foca = 0.4;
pos_fc   = 1.6;

Violin({nrmse_foca}, pos_foca, ...
    'ViolinColor', {C}, ...
    'HalfViolin', 'full', ...
    'ShowData', false, ...
    'ShowBox', true, ...
    'ViolinAlpha',{1}, ...
    'ShowMedian', true);

Violin({nrmse_fc}, pos_fc, ...
    'ViolinColor', {C}, ...
    'HalfViolin', 'full', ...
    'ShowData', false, ...
    'ShowBox', true, ...
    'ViolinAlpha',{1}, ...
    'ShowMedian', true);

set(gca, 'Box','off', 'LineWidth',1, ...
         'TickDir','out', 'TickLength',[.005 .005], ...
         'FontName','Arial','FontSize',12, ...
         'XTick',[pos_foca, pos_fc], ...
         'XTickLabel',{'FOCA','FC'}, ...
         'YTick',0:0.5:1.5, ...
         'YLim',[-0.1 1.4]);
ylabel('Normalized RMSE');
title('Stability Comparison (GSR vs. non-GSR)');
axis square;

%%
print('./results/StabilityComparison/FOCA_FC_Stability_Violin.png','-r300','-dpng');
disp('--- Stability comparison complete ---');