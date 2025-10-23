%% ================================================================
% Script: NeonatalAgePrediction_SVR_FOCA.m
% Purpose:
%   Validate predictive validity of FOCA (Functional Orthogonal 
%   Connectivity Architecture) features in predicting chronological 
%   age during the neonatal period using linear SVR.
%
% Description:
%   - FOCA features (GSR-preprocessed, individual-level dHCP data) 
%     are used as predictors.
%   - Sex, mean FD, and scan–birth intervals are regressed out 
%     from both FOCA features and target (chronological age).
%   - A 10-fold cross-validation linear SVR model is trained to 
%     predict postmenstrual age at scan.
%   - Prediction accuracy is defined as the Pearson correlation 
%     between predicted and actual age in held-out folds.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;
addpath('./dependencies');

%% ---------------------- Path Settings ----------------------------
PredictionFolder = './results/AgePrediction';
FeatureFolder = './data/dHCP/FOCA/';
AgeInfoFile = './data/dHCP/PredictionAgeScore.mat';
GroupMaskFile = './data/dHCP/Group_FC_dHCP_GSR.mat';

FoldQuantity = 10;
Pre_Method = 'Normalize';
C_Parameter = 1;
Permutation_Flag = 0;

%% ---------------------- Load FOCA Features -----------------------
disp('--- Loading FOCA (individual-level, dHCP GSR) features ---');
load([FeatureFolder 'FOCA_Individual_dHCP_GSR.mat']); 
FOCA_Individual = FOCA_Individual_dHCP_GSR;

% Remove subjects with zero maps
subject_sum = sum(FOCA_Individual, 2);
valid_idx = find(subject_sum ~= 0);
FOCA_Individual = FOCA_Individual(valid_idx, :);

% Replace any placeholder values (e.g., 1) with 0
FOCA_Individual(FOCA_Individual == 1) = 0;
clear subject_sum FOCA_Individual_dHCP_GSR

%% ---------------------- Load Age and Covariates ------------------
disp('--- Loading age and covariates ---');
load(AgeInfoFile, 'PredictionAgeScore');

Age = PredictionAgeScore.scan_age(valid_idx);
Sex = double(cellfun(@(x) x, PredictionAgeScore.sex(valid_idx, 2)));
MeanFD = PredictionAgeScore.MeanFD(valid_idx);
ScanBirthInterval = PredictionAgeScore.scan_age(valid_idx) - PredictionAgeScore.birth_age(valid_idx);
Covariates = [Sex, MeanFD, ScanBirthInterval];

%% ---------------------- Regress Out Covariates -------------------
disp('--- Regressing covariates from age and FOCA features ---');
% Remove covariate effects from age
stats = regstats(Age, Covariates, 'linear');
Age_resid = stats.beta(1) + stats.r;

% Load group-level FC to mask positive (retain negative) connections
load(GroupMaskFile, 'uCi_Group_FC');
neg_mask = (uCi_Group_FC <= 0);
FOCA_Individual(:, ~neg_mask(:)') = 0;

% Regress covariates from each FOCA feature
for j = 1:size(FOCA_Individual, 2)
    [b,~,~,~,~] = regress(FOCA_Individual(:, j), [ones(size(Covariates,1),1), Covariates]);
    FOCA_Individual(:, j) = FOCA_Individual(:, j) - Covariates*b(2:end);
end

clear b stats j

%% ---------------------- Run SVR Prediction -----------------------
disp('--- Running 10-fold SVR prediction (FOCA) ---');
ResultantFolder = [PredictionFolder '/SVR_Age_FOCA_dHCP_GSR_10CV/'];
if ~exist(ResultantFolder, 'dir')
    mkdir(ResultantFolder);
end

Prediction = SVR_NFolds_Sort( ...
    FOCA_Individual, Age_resid, FoldQuantity, ...
    Pre_Method, C_Parameter, Permutation_Flag, ResultantFolder);

fprintf('Prediction Results:\n');
fprintf('  r = %.3f\n', Prediction.r_value);
fprintf('  p = %.3e\n', Prediction.p_value);

% Save prediction results
save([ResultantFolder 'Prediction_FOCA_dHCP_GSR.mat'], ...
    'Prediction', 'FOCA_Individual', 'Age_resid');

%% ---------------------- Display Summary --------------------------
disp('--- FOCA (dHCP GSR) age prediction complete ---');
disp(['r = ' num2str(Prediction.r_value) ', p = ' num2str(Prediction.p_value)]);
