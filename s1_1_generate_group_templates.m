%% ================================================================
% Script: generate_group_templates.m
% Purpose:
%   Generate group-level functional network templates from individual 
%   seed-based correlation maps using the MSC dataset.
%
% Description:
%   - Loads a 20-network reference parcellation
%   - Computes seed-based correlation maps for each network per subject
%   - Averages Fisher-transformed correlations across subjects
%   - Outputs group-level MSC-derived templates
%
% Note:
%   The resulting individualized networks will be thresholded at z > 1 
%   (corresponding to the top ~15.9% of connections) during the 
%   individualized network generation step.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;

%% --- Dependencies ---
addpath(genpath('./dependence/cifti-matlab-master_kraus'));
addpath(genpath('./dependence/compare_matrices_to_assign_networks-main'));

%% --- Paths & Parameters ---
OutDir = './templates/';
SubResDir = './subjects/';
if ~exist(OutDir, 'dir'); mkdir(OutDir); end

% Subject list (MSC dataset)
Subject_all = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};

% Reference network definition (20 networks)
RefFile = fullfile('./templates/', 'Prior_Lynch20_Network.dtseries.nii');
consen = ft_read_cifti_mod(RefFile);
consen.data = consen.data(1:59412);  % cortical surface vertices only

network_names = { ...
    'Default_Parietal','Default_Anterolateral','Default_Dorsolateral','Default_Retrosplenial',...
    'Visual_Lateral','VentralStream','Visual_V5','Visual_V1',...
    'Frontoparietal','DorsalAttention','Premotor/DorsalAttentionII','Language',...
    'Salience','CinguloOpercular/Action-mode','MedialParietal',...
    'Somatomotor_Hand','Somatomotor_Face','Somatomotor_Foot',...
    'Auditory','SomatoCognitiveAction'};

%% --- Step 1. Compute seed-based correlation maps per subject ---
disp('Computing subject-level seed-based correlation maps...');

for i = 1:length(Subject_all)
    Subject = Subject_all{i};
    PfmDir = fullfile(SubResDir, 'pfm', Subject);

    % Load subject BOLD data
    boldFile = fullfile(PfmDir, ['sub-' Subject '_task-rest_concatenated_smoothed2.55_32k_fsLR.dtseries.nii']);
    TR = ft_read_cifti_mod(boldFile);

    % Compute correlation for each network
    corrs = nan(59412, length(network_names));
    for j = 1:length(network_names)
        disp(['  Processing network ' num2str(j) ': ' network_names{j}]);
        inds = (consen.data == j);
        subNetAvg = nanmean(TR.data(inds,:), 1);

        for voxel = 1:size(TR.data, 1)
            goodvox = ~isnan(TR.data(voxel, :));
            corrs(voxel, j) = paircorr_mod(subNetAvg(goodvox)', TR.data(voxel, goodvox)')';
        end
    end

    seedmapsTR{i} = corrs;
    clear TR corrs subNetAvg
end

%% --- Step 2. Average across subjects (Fisher z-transform) ---
disp('Averaging across subjects...');
subs = length(Subject_all);
nVert = 59412;
nNet = length(network_names);
Prior = zeros(nVert, nNet);

for j = 1:nNet
    grpNetAve = zeros(nVert, 1);
    for i = 1:subs
        grpNetAve = grpNetAve + atanh(seedmapsTR{i}(1:nVert, j)); % Fisher transform
    end
    grpNetAve = grpNetAve ./ subs;
    Prior(:, j) = grpNetAve;
end

%% --- Step 3. Save outputs ---
disp('Saving MSC-derived group templates...');
temp = ft_read_cifti_mod(boldFile); % use header from a subject file
temp.data = zeros(size(temp.data,1), nNet);
temp.data(1:nVert, :) = Prior;

OutTemplate = fullfile(OutDir, 'Prior_MSC_Template_Lynch20.dtseries.nii');
ft_write_cifti_mod(OutTemplate, temp);
save(fullfile(OutDir, 'Prior_MSC_Lynch20_fisher.mat'), 'Prior');

disp('MSC-derived group-level templates generated successfully!');
