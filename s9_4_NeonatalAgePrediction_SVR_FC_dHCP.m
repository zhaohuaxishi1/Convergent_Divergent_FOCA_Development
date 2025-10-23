%% ================================================================
% Script: NeonatalAgePrediction_SVR_FC_dHCP_GSR.m
% Purpose:
%   Evaluate the predictive validity of conventional functional 
%   connectivity (FC) features in predicting neonatal chronological age 
%   using linear SVR with 10-fold cross-validation.
%
% Description:
%   - Individual-level FC features (GSR-preprocessed, dHCP dataset) 
%     are used as predictors.
%   - Sex, mean framewise displacement (FD), and scan–birth intervals 
%     are regressed out from both FC features and the target variable (age).
%   - A linear SVR model with 10-fold CV is applied.
%   - Prediction accuracy is defined as the Pearson correlation between 
%     predicted and actual age across held-out folds.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;
addpath('./dependencies/SVR');

%% ---------------------- Path Settings ----------------------------
PredictionFolder = './results/AgePrediction';
FeatureFolder = './data/dHCP/FC/';
AgeInfoFile = './data/dHCP/PredictionAgeScore.mat';
GroupMaskFile = './data/dHCP/Group_FC_dHCP_GSR.mat';

FoldQuantity = 10;
Pre_Method = 'Normalize';
C_Parameter = 1;
Permutation_Flag = 0;

%% ---------------------- Load FC Features -------------------------
disp('--- Loading individual-level FC (dHCP GSR) features ---');
load([FeatureFolder 'FC_Individual_dHCP_GSR.mat']); 
FC_Individual = FC_Individual_dHCP_GSR;

% Remove subjects with zero connectivity
subject_sum = sum(FC_Individual, 2);
valid_idx = find(subject_sum ~= 0);
FC_Individual = FC_Individual(valid_idx, :);

% Replace placeholder values (if any)
FC_Individual(FC_Individual == 1) = 0;
clear subject_sum FC_Individual_dHCP_GSR

%% ---------------------- Load Age and Covariates ------------------
disp('--- Loading age and covariates ---');
load(AgeInfoFile, 'PredictionAgeScore');

Age = PredictionAgeScore.scan_age(valid_idx);
Sex = double(cellfun(@(x) x, PredictionAgeScore.sex(valid_idx, 2)));
MeanFD = PredictionAgeScore.MeanFD(valid_idx);
ScanBirthInterval = PredictionAgeScore.scan_age(valid_idx) - PredictionAgeScore.birth_age(valid_idx);
Covariates = [Sex, MeanFD, ScanBirthInterval];

%% ---------------------- Regress Out Covariates -------------------
disp('--- Regressing covariates from age and FC features ---');
% Regress covariate effects from age
stats = regstats(Age, Covariates, 'linear');
Age_resid = stats.beta(1) + stats.r;

% Load group-level FC matrix to identify positive connections
load(GroupMaskFile, 'uCi_Group_FC');
pos_mask = (uCi_Group_FC > 0);
FC_Individual(:, pos_mask(:)') = 0; % Retain negative FC only

% Regress covariates from each FC feature
for j = 1:size(FC_Individual, 2)
    [b,~,~,~,~] = regress(FC_Individual(:, j), [ones(size(Covariates,1),1), Covariates]);
    FC_Individual(:, j) = FC_Individual(:, j) - Covariates*b(2:end);
end

clear b stats j

%% ---------------------- Run SVR Prediction -----------------------
disp('--- Running 10-fold SVR prediction (FC) ---');
ResultantFolder = [PredictionFolder '/SVR_Age_FC_dHCP_GSR_10CV/'];
if ~exist(ResultantFolder, 'dir')
    mkdir(ResultantFolder);
end

Prediction = SVR_NFolds_Sort( ...
    FC_Individual, Age_resid, FoldQuantity, ...
    Pre_Method, C_Parameter, Permutation_Flag, ResultantFolder);

fprintf('Prediction Results:\n');
fprintf('  r = %.3f\n', Prediction.r_value);
fprintf('  p = %.3e\n', Prediction.p_value);

% Save prediction results
save([ResultantFolder 'Prediction_FC_dHCP_GSR.mat'], ...
    'Prediction', 'FC_Individual', 'Age_resid');

%% ---------------------- Display Summary --------------------------
disp('--- FC (dHCP GSR) age prediction complete ---');
disp(['r = ' num2str(Prediction.r_value) ', p = ' num2str(Prediction.p_value)]);
