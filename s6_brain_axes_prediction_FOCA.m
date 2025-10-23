%% ================================================================ 
% Script: brain_axes_prediction_FOCA.m
% Purpose:
%   Test whether internetwork covariance (FOCA) is shaped by 
%   fundamental cortical organizational features using multiple 
%   linear regression with 10-fold cross-validation.
%
% Description:
%   - Incorporate nine cortical neurobiological axes from neuromaps.
%   - Derive network-level feature alignment matrices.
%   - Predict FOCA matrix via MultipleRegression_NFolds_Sort function.
%   - Estimate relative contribution of each cortical feature.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
%
% ================================================================

clc; clear;
addpath(genpath('./dependence/customcolormap'));
addpath(genpath('./dependence/cifti-matlab-master_kraus'));
addpath(genpath('./dependence/MultipleRegression_NFolds_Sort')); % ensure function path is added

%% --- Step 1. Load cortical neuromaps (nine cortical features) ---
disp('--- Loading cortical neurobiological axes ---');
neuromap_dir = './data/S-A_ArchetypalAxis-main/FSLRVertex/';
BrainMap = readmatrix([neuromap_dir 'brainmaps_fslr.csv'], 'NumHeaderLines', 1); 
BrainMap(isnan(BrainMap)) = 0;

% Use the sensorimotor–association (S-A) axis for orientation consistency
SA_axis = ft_read_cifti_mod([neuromap_dir 'SensorimotorAssociation_Axis.dscalar.nii']);
SA_axis = SA_axis.data;

flip_idx = [];
for i = 1:size(BrainMap,2)
    [r,p] = corr(SA_axis, BrainMap(:,i));
    if r < 0, flip_idx = [flip_idx i]; end
end
BrainMap(:,flip_idx) = -BrainMap(:,flip_idx);
clear SA_axis flip_idx r p

%% --- Step 2. Load group-level individualized topographies (MSC dataset) ---
disp('--- Loading group-level individualized topography ---');
load('./data/Individual_Topography_Lynch20_MSC_all.mat');
nNet = 20;

uCi_Group_Topography = zeros(size(uCi_individual_network_all{1}));
for i = 1:length(uCi_individual_network_all)
    uCi_Group_Topography = uCi_Group_Topography + uCi_individual_network_all{i};
end
uCi_Group_Topography = uCi_Group_Topography ./ length(uCi_individual_network_all);
[~, Group_Parcellation] = max(uCi_Group_Topography,[],2);
uCi = unique(nonzeros(Group_Parcellation));

%% --- Step 3. Compute network-level alignment matrices for each cortical axis ---
disp('--- Computing network-level alignment matrices ---');
num_axes = size(BrainMap, 2);
Lynch20_Map_Matrix_all = cell(1, num_axes);

for k = 1:num_axes
    feature_values = zeros(length(uCi), 1);
    for i = 1:length(uCi)
        feature_values(i) = mean(BrainMap(Group_Parcellation==uCi(i), k), 'omitnan');
    end
    feature_values = rescale(feature_values, -1, 1);
    
    feature_matrix = zeros(length(uCi));
    for i = 1:length(uCi)
        for j = 1:length(uCi)
            if i == j
                feature_matrix(i,j) = 1;
            else
                feature_matrix(i,j) = 1 - abs(feature_values(i) - feature_values(j));
            end
        end
    end
    Lynch20_Map_Matrix_all{k} = feature_matrix;
end
save('./results/FOCA/Lynch20_Map_Matrix_all_neuromaps.mat', 'Lynch20_Map_Matrix_all');

%% --- Step 4. Load FOCA matrix ---
disp('--- Loading FOCA matrix ---');
load('./results/FOCA/MSC_FOCA_GroupLevel_Matrix.mat', 'corrMatrix');
Y = corrMatrix(:);

%% --- Step 5. Prepare predictors for regression ---
disp('--- Preparing cortical feature predictors ---');
% Order of nine cortical features:
% AH, EH, AS, AG, CBF, GE, NS, EX, CT
X = [Lynch20_Map_Matrix_all{1}(:),Lynch20_Map_Matrix_all{3}(:), ...
     Lynch20_Map_Matrix_all{4}(:),Lynch20_Map_Matrix_all{5}(:), ...
     Lynch20_Map_Matrix_all{6}(:),Lynch20_Map_Matrix_all{7}(:), ...
     Lynch20_Map_Matrix_all{8}(:),Lynch20_Map_Matrix_all{9}(:), ...
     Lynch20_Map_Matrix_all{10}(:)];

%% --- Step 6. Run regression using MultipleRegression_NFolds_Sort ---
disp('--- Running multiple regression (10-fold CV) ---');
PredictionFolder = './results/FOCA/BrainAxes_Prediction/';
if ~exist(PredictionFolder, 'dir'), mkdir(PredictionFolder); end

Neuromaps = X;
Inter_network_matrix = Y;
FoldQuantity = 10;
Pre_Method = 'Normalize';
C_Parameter = 1; 
Permutation_Flag = 0; % no permutation test
ResultantFolder = [PredictionFolder 'MultipleRegression_10CV_Sort/'];
if ~exist(ResultantFolder, 'dir'), mkdir(ResultantFolder); end

Prediction = MultipleRegression_NFolds_Sort( ...
    Neuromaps, Inter_network_matrix, FoldQuantity, ...
    Pre_Method, Permutation_Flag, ResultantFolder);

% Save prediction results
save([ResultantFolder 'FOCA_BrainAxes_Prediction_Result.mat'], ...
    'Prediction', 'Neuromaps', 'Inter_network_matrix');

%% --- Step 7. Compute and visualize relative contribution ---
disp('--- Computing contribution weights ---');
X = Neuromaps;
Y = Inter_network_matrix;
stat = regstats(Y,X);
group_pre_FC_r2 = stat.adjrsquare;
% disp(stat.beta);

loadings = stat.beta(2:end);
contributionPercent = 100 * (loadings.^2) / sum(loadings.^2);
for i = 1:9
    contributionPercent(i,2) = i;
end
disp(contributionPercent);
[V_median_sorted,V_median_sorted_index] = sort(contributionPercent(:,1), 'descend');
V_median_result = contributionPercent(V_median_sorted_index,:);
dataset = V_median_result(:,1);

%% --- Step 8. Visualization (standardized publication style) ---
disp('--- Visualizing cortical feature contributions ---');

% Define colors
C1 = [83,20,112]/255;  % purple
C2 = [237,177,32]/255; % gold

% Figure settings
figureUnits = 'centimeters';
figureWidth = 12;
figureHeight = 10;
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);
hold on;

% Assign bar colors
colors = [C1; C1; C2; C2; C2; C2; C2; C2; C2];
for i = 1:9
    GO = bar(i, dataset(i), 0.75);
    set(GO(1),'facecolor',colors(i,:));
end

% Axis style
set(gca, 'Box', 'off', 'LineWidth', 1, ...
    'XGrid', 'off', 'YGrid', 'off', ...
    'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
    'YTick', 0:20:60, 'Ylim', [0 60], ...
    'XTick', 1:1:9, 'Xlim', [0 10], ...
    'Xticklabel', {'AG' 'EX' 'GE' 'NS' 'AS' 'EH' 'CBF' 'AH' 'CT'}, ...
    'Yticklabel', [0:20:60]);

hYLabel = ylabel('Relative contribution (%)');
set(gca,'fontname','arial','fontsize',12);
set(gca,'Color',[1 1 1]);

% Export figure
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figureWidth figureHeight]);
fileout = [ResultantFolder 'FOCA_BrainAxes_Contribution'];
print(figureHandle,[fileout,'.png'],'-r300','-dpng');

disp('--- Brain axes regression and contribution visualization complete ---');
