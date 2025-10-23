%% ================================================================
% Script: compute_foca_matrix.m
% Purpose:
%   Construct Functional Topography Covariance (FOCA) matrices 
%   based on individualized functional topographies.
%
% Description:
%   - Loads subject-level individualized topography maps
%   - Computes Pearson correlations between spatial maps of all network pairs
%   - Produces one FOCA matrix per subject
%   - Saves FOCA matrices and optional heatmap visualizations
%
% Interpretation:
%   Positive correlations reflect cooperative interactions,
%   whereas negative correlations indicate competitive interactions.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;

%% --- Dependencies ---
addpath(genpath('./dependence/customcolormap'));
addpath(genpath('./dependence/cifti-matlab-master'));

%% --- Paths & Parameters ---
Subject_all = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
DataDir = './data/';                 % folder containing Individual_Topography_Lynch20_MSC_all.mat
OutDir  = './results/FOCA/';            % output directory for FOCA matrices
if ~exist(OutDir, 'dir'); mkdir(OutDir); end

disp('Computing FOCA matrices for all subjects...');

%% --- Step 1. Load subject-level individualized topographies ---
load(fullfile(DataDir, 'Individual_Topography_Lynch20_MSC_all.mat'));  % variable: uCi_individual_network_all
nSubjects = length(Subject_all);

% Define network labels (20 networks)
network_labels = { ...
    'DMN-Par', 'DMN-Ant', 'DMN-Dor', 'DMN-Ret', ...
    'Vis-Lat', 'Vis-Str', 'Vis-V5', 'Vis-V1', ...
    'FP', 'DAN', 'PreMot', 'Lang', ...
    'Sal', 'CO/AMN', 'PMN', ...
    'Hand', 'Face', 'Foot', 'Aud', 'SCAN'};

%% --- Step 2. Compute FOCA (network-by-network correlation) ---
for i = 1:nSubjects
    Subject = Subject_all{i};
    disp(['Processing subject: ' Subject]);

    if isempty(uCi_individual_network_all{i})
        warning(['No data found for subject ' Subject ', skipping...']);
        continue;
    end

    data = uCi_individual_network_all{i};  % size: vertices × networks

    % Compute Pearson correlation between networks (20 × 20)
    FOCA_matrix = corr(data);

    % --- Save FOCA matrix ---
    save(fullfile(OutDir, ['FOCA_matrix_' Subject '.mat']), 'FOCA_matrix');

    %% --- Optional: Visualization ---
    figure('Name', ['FOCA - ' Subject], 'NumberTitle', 'off', ...
           'Units', 'inches', 'Position', [0 0 9 8]);
    imagesc(FOCA_matrix);
    colormap(customcolormap_preset('red-white-blue'));
    caxis([-0.6 1]);
    colorbar;
    cb = colorbar;
    cb.Ticks = -0.6:0.4:1;

    ax = gca;
    set(ax, 'XTick', 1:20, 'YTick', 1:20, ...
        'XTickLabel', network_labels, 'YTickLabel', network_labels, ...
        'TickLabelInterpreter', 'none', 'FontName', 'Arial', 'FontSize', 10);
    xtickangle(90);
    axis square;

    % Save figure
    print(gcf, fullfile(OutDir, ['FOCA_matrix_' Subject '.png']), '-r300', '-dpng');
    close(gcf);
end

disp('All FOCA matrices computed successfully!');
