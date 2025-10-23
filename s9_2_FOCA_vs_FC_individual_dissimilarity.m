%% ================================================================
% Script: FOCA_vs_FC_individual_dissimilarity.m
% Purpose:
%   Evaluate whether FOCA provides individual-specific information
%   beyond conventional FC by comparing subject-level dissimilarity
%   (1 - Pearson’s r) between group- and individual-level maps of
%   negative internetwork connections.
%
% Description:
%   - Compute group-level FOCA/FC matrices (GSR).
%   - Compute individual-level FOCA/FC matrices.
%   - Restrict to negative group-level connections.
%   - For each subject: dissimilarity = 1 - corr(group, individual).
%   - Compare FOCA vs FC dissimilarity using paired t-test.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;
addpath('./utilities');
addpath('./dependencies');
addpath('./visualization');

%% --- Step 1. Load FOCA (Functional Topographic Covariance) data ---
disp('--- Loading FOCA data ---');
load('./data/HCP/FOCA_Group_HCP_GSR.mat');              % group-level FOCA
FOCA_Group_GSR = corr(FOCA_Group);                      % convert to correlation (20×20)

load('./data/HCP/FOCA_Individual_HCP_GSR.mat'); % individual-level FOCA
FOCA_Individual_GSR = FOCA_Individual_HCP_GSR;
clear FOCA_Group FOCA_Individual_HCP_GSR

% Remove zero subjects
valid_idx = find(sum(FOCA_Individual_GSR,2)~=0);
FOCA_Individual_GSR = FOCA_Individual_GSR(valid_idx,:);

% Identify negative connections only
neg_mask = FOCA_Group_GSR < 0;
neg_idx  = find(neg_mask);

% Compute dissimilarity for each subject
nSub = size(FOCA_Individual_GSR,1);
FOCA_r = nan(nSub,1);

for s = 1:nSub
    indi_mat = reshape(FOCA_Individual_GSR(s,:), 20, 20);
    [FOCA_r(s), ~] = corr(FOCA_Group_GSR(neg_idx), indi_mat(neg_idx), 'rows', 'complete');
end

FOCA_dissimilarity = 1 - FOCA_r;  % subject-level dissimilarity

fprintf('FOCA mean dissimilarity = %.3f ± %.3f\n', mean(FOCA_dissimilarity), std(FOCA_dissimilarity));

%% --- Step 2. Load conventional FC data ---
disp('--- Loading FC data ---');
load('./data/HCP/FC_Group_HCP_GSR.mat');             % group-level FC
FC_Group_GSR = FC_Group;

load('./data/HCP/FC_Individual_HCP_GSR.mat');        % individual-level FC
FC_Individual_GSR = FC_Individual_HCP_GSR;
clear FC_Group FC_Individual_HCP_GSR

% Remove zero subjects
valid_idx = find(sum(FC_Individual_GSR,2)~=0);
FC_Individual_GSR = FC_Individual_GSR(valid_idx,:);

% Compute dissimilarity for FC
nSub = size(FC_Individual_GSR,1);
FC_r = nan(nSub,1);

for s = 1:nSub
    indi_mat = reshape(FC_Individual_GSR(s,:), 20, 20);
    [FC_r(s), ~] = corr(FC_Group_GSR(neg_idx), indi_mat(neg_idx), 'rows', 'complete');
end

FC_dissimilarity = 1 - FC_r;

fprintf('FC mean dissimilarity = %.3f ± %.3f\n', mean(FC_dissimilarity), std(FC_dissimilarity));

%% --- Step 3. Statistical comparison (paired t-test) ---
disp('--- Comparing FOCA vs FC individual dissimilarity ---');
[~, p_value, ~, stats] = ttest2(FOCA_dissimilarity, FC_dissimilarity);
T_value = stats.tstat;
Cohen_d_value = computeCohen_d(FOCA_dissimilarity, FC_dissimilarity);

fprintf('Paired t-test: T = %.3f, p = %.3e, Cohen''s d = %.3f\n', ...
        T_value, p_value, Cohen_d_value);

save('./results/IndividualDissimilarity/FOCA_FC_Dissimilarity_Stats.mat', ...
     'FOCA_dissimilarity', 'FC_dissimilarity', 'T_value', 'p_value', 'Cohen_d_value');

%% --- Step 4. Visualization ---
disp('--- Visualizing FOCA vs FC dissimilarity ---');
load('./data/Priors/priors.mat');
C = Priors.NetworkColors;

figure('Units','centimeters','Position',[0 0 6 7]);
hold on;

pos_foca = 0.4;
pos_fc   = 1.6;

Violin({FOCA_dissimilarity}, pos_foca, ...
    'ViolinColor', {C(13,:)}, ...
    'HalfViolin', 'full', ...
    'ShowData', false, ...
    'ShowBox', true, ...
    'ViolinAlpha',{1}, ...
    'ShowMedian', true);

Violin({FC_dissimilarity}, pos_fc, ...
    'ViolinColor', {C(13,:)}, ...
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
         'YTick',0:0.2:0.8, ...
         'YLim',[0 0.9]);
ylabel('Individual dissimilarity (1 - r)');
title('FOCA vs FC negative-connection dissimilarity');
axis square;

print('./results/IndividualDissimilarity/FOCA_FC_Dissimilarity_Violin.png','-r300','-dpng');
disp('--- Individual dissimilarity comparison complete ---');
