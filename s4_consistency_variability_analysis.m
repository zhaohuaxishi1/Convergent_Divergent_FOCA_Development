%% ================================================================
% Script: consistency_variability_analysis.m
% Purpose:
%   Assess cross-subject consistency (MSC dataset) and inter-individual 
%   variability (HCP dataset) of FOCA matrices.
%
% Description:
%   - Part I: MSC dataset — compute pairwise correlations between 
%     participants’ FOCA matrices to assess cross-subject consistency.
%   - Part II: HCP dataset — compute Median Absolute Deviation (MAD) 
%     across participants as a robust index of inter-individual variability.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
%
% Notes:
%   • Cross-subject consistency (MSC): mean pairwise r of FOCA matrices.
%   • Inter-individual variability (HCP): MAD across subjects (robust to outliers).
%   • Because of its large sample size, HCP data provide a more reliable 
%     estimation of inter-individual variability compared to MSC.
% ================================================================

clc; clear;

%% ================================================================
% PART I — MSC: Cross-subject FOCA consistency
% ================================================================

addpath(genpath('./dependence/customcolormap'));

DataDir_MSC = './results/FOCA/';
OutDir_MSC  = './results/FOCA_Analysis_MSC/';
if ~exist(OutDir_MSC, 'dir'); mkdir(OutDir_MSC); end

Subject_all = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
nSub = length(Subject_all);
nNet = 20;

disp('--- [MSC] Assessing cross-subject consistency of FOCA matrices ---');

%% --- Load FOCA matrices for MSC subjects ---
FOCA_all = zeros(nNet, nNet, nSub);
for i = 1:nSub
    file = fullfile(DataDir_MSC, ['FOCA_matrix_' Subject_all{i} '.mat']);
    if isfile(file)
        tmp = load(file, 'FOCA_matrix');
        FOCA_all(:, :, i) = tmp.FOCA_matrix;
    else
        warning(['Missing FOCA matrix for subject ' Subject_all{i}]);
        FOCA_all(:, :, i) = nan;
    end
end

%% --- Compute pairwise correlation (cross-subject consistency) ---
Sub_r = nan(nSub, nSub);
for i = 1:nSub
    for j = 1:nSub
        if any(isnan(FOCA_all(:, :, i)), 'all') || any(isnan(FOCA_all(:, :, j)), 'all')
            continue;
        end
        mat1 = FOCA_all(:, :, i);
        mat2 = FOCA_all(:, :, j);
        Sub_r(i, j) = corr(mat1(:), mat2(:), 'rows', 'pairwise');
    end
end

mean_r = mean(Sub_r(~isnan(Sub_r) & triu(true(size(Sub_r)), 1)), 'all');
std_r  = std(Sub_r(~isnan(Sub_r) & triu(true(size(Sub_r)), 1)), 0, 'all');
fprintf('\n[MSC] Cross-subject FOCA consistency:\nMean r = %.2f, SD = %.2f\n', mean_r, std_r);

save(fullfile(OutDir_MSC, 'FOCA_consistency_matrix_MSC.mat'), 'Sub_r', 'mean_r', 'std_r');

%% --- Visualize MSC cross-subject consistency ---
figure('Name','MSC FOCA Cross-subject Consistency','NumberTitle','off','Units','inches','Position',[0 0 8 7]);
imagesc(Sub_r);
colormap(jet);
colorbar;
caxis([0 1]);
cb = colorbar; cb.Ticks = 0:0.2:1;
set(gca,'XTick',1:nSub,'YTick',1:nSub,'XTickLabel',Subject_all,'YTickLabel',Subject_all,...
    'TickLabelInterpreter','none','FontSize',10,'FontName','Arial');
xtickangle(90); axis square;
title('MSC FOCA Cross-subject Consistency');
print(gcf, fullfile(OutDir_MSC, 'FOCA_consistency_MSC_heatmap.png'), '-r300', '-dpng');
close(gcf);

%% ================================================================
% PART II — HCP: Inter-individual FOCA variability (MAD)
% ================================================================

disp('--- [HCP] Computing inter-individual FOCA variability (MAD) ---');

DataFile_HCP = './data/HCP/FOCA_Individual_HCP_GSR.mat';
OutDir_HCP   = './results/FOCA_Analysis_HCP/';
if ~exist(OutDir_HCP, 'dir'); mkdir(OutDir_HCP); end

load(DataFile_HCP);   % variable: HCP_FOCA_AllSubjects
SubjectsData = FOCA_Individual_HCP_GSR;          % Each row = subject, each column = flattened FOCA vector

% Clean and filter
SubjectsData(SubjectsData == 1) = 0;  % remove spurious values
AtlasLoading_Sum = sum(SubjectsData, 2);
SubjectsData = SubjectsData(AtlasLoading_Sum > 0, :);

%% --- Compute MAD across subjects ---
madRow = mad(SubjectsData, 1, 1);  % median absolute deviation per feature
madMat = reshape(madRow, nNet, nNet);

meanMAD = mean(madMat(:), 'omitnan');
stdMAD  = std(madMat(:), 'omitnan');
fprintf('\n[HCP] FOCA inter-individual variability:\nMean MAD = %.3f, SD = %.3f\n', meanMAD, stdMAD);

save(fullfile(OutDir_HCP, 'HCP_FOCA_variability_MAD.mat'), 'madMat', 'meanMAD', 'stdMAD');

%% --- Compute intra- and inter-network variability ---
corrMatrixLabel = [1,1,1,1,2,2,2,2,3,4,4,5,6,7,8,9,9,9,10,11];
uCi = unique(nonzeros(corrMatrixLabel));
uCi_FC = zeros(length(uCi));

for i = 1:length(uCi)
    for j = 1:length(uCi)
        uCi_FC(i,j) = nanmean(nanmean(madMat(corrMatrixLabel==uCi(i), corrMatrixLabel==uCi(j))));
    end
end

diag_mean     = mean(diag(uCi_FC));           % intra-network variability
nondiag_mean  = mean(uCi_FC(~eye(size(uCi_FC)))); % inter-network variability

fprintf('\n[HCP] Intra-network MAD = %.3f | Inter-network MAD = %.3f\n', diag_mean, nondiag_mean);

save(fullfile(OutDir_HCP, 'HCP_FOCA_variability_summary.mat'), ...
     'uCi_FC', 'diag_mean', 'nondiag_mean', 'corrMatrixLabel');

%% --- Visualize HCP MAD matrix ---
network_labels = { ...
    'DMN Par','DMN Ant','DMN Dor','DMN Ret', ...
    'Vis Lat','Vis Str','Vis V5','Vis V1', ...
    'FP','DAN','PreMot','Lang','Sal','CO/AMN', ...
    'PMN','Hand','Face','Foot','Aud','SCAN'};

figure('Name','HCP FOCA Inter-individual Variability (MAD)', ...
       'NumberTitle','off','Units','inches','Position',[0 0 9 8]);
imagesc(madMat);
colormap(jet);
caxis([0 0.4]);
cb = colorbar; cb.Ticks = 0:0.1:0.4;
set(gca,'XTick',1:nNet,'YTick',1:nNet,...
    'XTickLabel',network_labels,'YTickLabel',network_labels,...
    'TickLabelInterpreter','none','FontName','Arial','FontSize',10);
xtickangle(90); axis square;
title('HCP FOCA Inter-individual Variability (MAD)');
print(gcf, fullfile(OutDir_HCP, 'HCP_FOCA_variability_MAD_heatmap.png'), '-r300', '-dpng');
close(gcf);

disp('--- FOCA consistency and variability analysis completed successfully ---');
