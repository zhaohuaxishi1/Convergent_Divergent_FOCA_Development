%% ================================================================
% Script: generate_individual_networks.m
% Purpose:
%   Generate individualized functional networks for each subject 
%   using the template-matching (TM) approach.
%
% Description:
%   - Applies MSC-derived group templates to each subject (MSC/HCP dataset)
%   - Computes whole-brain functional connectivity profiles per grayordinate
%   - Thresholds connectivity maps (z > 1; top ~15.9% of connections)
%   - Calculates η² values between each grayordinate’s connectivity 
%     profile and each network template
%   - Assigns each grayordinate to the network with the highest η² value
%     to generate individualized network topography (Supplementary Fig. S4)
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
%
% Note:
%   The same thresholding criterion (z > 1) as the group-level template 
%   is applied for consistent network assignment across individuals.
% ================================================================

clc; clear;

%% --- Dependencies ---
addpath(genpath('./dependence/compare_matrices_to_assign_networks-main'));
addpath(genpath('./dependence/cifti-matlab-master_kraus'));

%% --- Paths & Parameters ---
Subject_all = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
template_path = './outputs/Prior_MSC_Lynch20_fisher.mat';    % group-level MSC template
SubResDir = './subjects/';                                   % subject directories
wb_command = '/usr/local/workbench/bin_linux64/wb_command';   % Connectome Workbench executable
disp('Generating individualized functional networks using Template Matching...');

%% --- Main loop over subjects ---
for i = 1:length(Subject_all)
    Subject = Subject_all{i};
    disp(['Processing subject: ' Subject]);

    PfmDir = fullfile(SubResDir, 'pfm', Subject);
    if ~exist(PfmDir, 'dir')
        warning(['Skipping ' Subject ' (missing data folder)']);
        continue;
    end

    % Input: subject fMRI data
    dconn_filename = fullfile(PfmDir, ...
        ['sub-' Subject '_task-rest_concatenated_smoothed2.55_32k_fsLR.dtseries.nii']);
    data_type = 'dtseries';

    % Output: individualized functional network assignment
    output_cifti_name = 'Individualized_Functional_Network_Assignment_Lynch20';
    cifti_output_folder = PfmDir;

    % TM parameters
    transform_data = 'Convert_r_to_Fisher';
    make_cifti_from_results = 1;
    allow_overlap = 0;
    overlap_method = "";
    surface_only = 1;
    already_surface_only = 1;

    % Load subject mid-thickness surface (for surface projection)
    Subdir = fullfile('./subjects/MSC/derivatives/surface_pipeline', ...
        ['sub-' Subject], 'fs_LR_Talairach', 'fsaverage_LR32k');
    MidthickSurfs{1} = fullfile(Subdir, [Subject '.L.midthickness.32k_fs_LR.surf.gii']);
    MidthickSurfs{2} = fullfile(Subdir, [Subject '.R.midthickness.32k_fs_LR.surf.gii']);

    % --- Run Template Matching ---
    template_matching_RH_Lynch20( ...
        dconn_filename, data_type, template_path, transform_data, ...
        output_cifti_name, cifti_output_folder, wb_command, ...
        make_cifti_from_results, allow_overlap, overlap_method, ...
        surface_only, already_surface_only, MidthickSurfs);

    disp(['Subject ' Subject ' completed.']);
end

disp('All subjects processed successfully!');
