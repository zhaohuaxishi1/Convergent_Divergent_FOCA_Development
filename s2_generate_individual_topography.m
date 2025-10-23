%% ================================================================
% Script: generate_individual_topography.m
% Purpose:
%   Derive participant-specific functional network topography 
%   from individualized functional networks.
%
% Description:
%   - Loads individualized functional network labels (from template matching)
%   - Computes the mean time series of each individualized network
%   - Correlates each network’s mean time series with all cortical vertices
%   - Generates vertex-wise topographic maps per subject (Fisher z-transformed)
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;

%% --- Dependencies ---
addpath(genpath('./dependence/cifti-matlab-master'));
addpath(genpath('./dependence/compare_matrices_to_assign_networks-main'));

%% --- Paths & Parameters ---
Subject_all = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
SubResDir = './subjects/';
OutDir = './data/';
if ~exist(OutDir, 'dir'); mkdir(OutDir); end

disp('Generating participant-specific topographic maps...');

tic;
uCi_individual_network_all = cell(length(Subject_all), 1);

%% --- Step 1. Compute network-to-vertex correlations per subject ---
for k = 1:length(Subject_all)
    Subject = Subject_all{k};
    disp(['Processing subject: ' Subject]);

    PfmDir = fullfile(SubResDir, 'pfm', Subject);

    % --- Load individualized network labels (from template matching) ---
    indvLabelFile = fullfile(PfmDir, 'Individualized_Functional_Network_Assignment_Lynch20.dlabel.nii');
    if ~isfile(indvLabelFile)
        warning(['Label file not found for subject ' Subject '. Skipping...']);
        continue;
    end
    Ic = ft_read_cifti_mod(indvLabelFile);

    % --- Load subject’s fMRI time series ---
    boldFile = fullfile(PfmDir, ['sub-' Subject '_task-rest_concatenated_smoothed2.55_32k_fsLR.dtseries.nii']);
    if ~isfile(boldFile)
        warning(['BOLD file missing for subject ' Subject '. Skipping...']);
        continue;
    end
    C = ft_read_cifti_mod(boldFile);
    nCorticalVertices = nnz(C.brainstructure == 1) + nnz(C.brainstructure == 2);

    % --- Identify unique individualized networks ---
    uCi = unique(nonzeros(Ic.data(1:nCorticalVertices, 1)));

    % --- Preallocate matrix for topographic maps ---
    uCi_individual_network = zeros(nCorticalVertices, length(uCi));

    % --- Compute mean time series per network and correlate with all vertices ---
    for j = 1:length(uCi)
        inds = (Ic.data(1:nCorticalVertices, 1) == uCi(j));
        subNetAvg = nanmean(C.data(inds, :), 1);
        for voxel = 1:nCorticalVertices
            goodvox = ~isnan(C.data(voxel, :));
            uCi_individual_network(voxel, j) = paircorr_mod(subNetAvg(goodvox)', C.data(voxel, goodvox)')';
        end
    end

    % --- Fisher r-to-z transformation ---
    uCi_individual_network(uCi_individual_network >= 1) = 0.9999;
    uCi_individual_network = atanh(uCi_individual_network);

    % --- Save subject-level results ---
    save(fullfile(PfmDir, 'Individual_Topography_Lynch20.mat'), 'uCi_individual_network');
    O = C;
    O.data = zeros(size(O.data, 1), size(uCi_individual_network, 2));
    O.data(1:nCorticalVertices, :) = uCi_individual_network;

    OutFile = fullfile(PfmDir, 'Individual_Topography_Lynch20.dtseries.nii');
    ft_write_cifti_mod(OutFile, O);

    uCi_individual_network_all{k} = uCi_individual_network;

    disp(['Subject ' Subject ' topographic maps generated.']);
    clear C O Ic uCi uCi_individual_network
end

%% --- Step 2. Save group collection ---
save(fullfile(OutDir, 'Individual_Topography_Lynch20_MSC_all.mat'), 'uCi_individual_network_all');

toc;
disp('All participant-specific topographic maps generated successfully.');
