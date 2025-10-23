%% ================================================================
% Script: behavior_prediction_neonatal_FOCA_GLM.m
% Purpose:
%   Test whether neonatal internetwork covariance (FOCA matrix)
%   predicts 18-month developmental outcomes (Bayley-III: Cognitive,
%   Language, Motor).
%
% Description:
%   - Each network’s FOCA profile (20×20) used as predictors.
%   - Dependent variables: Bayley-III composite scores (Cognitive, Language, Motor).
%   - Covariates: sex, mean FD, scan-birth interval, PMA at scan.
%   - Multiple linear regression applied for each network separately.
%   - Prediction accuracy evaluated via Pearson correlation between
%     actual and predicted behavioral scores.
%
% Reference:
%   Zhao et al., "Convergent and divergent functional topographies 
%   of individualized human brain network and their developmental origins"
% ================================================================

clc; clear;

%% --- Step 1. Load data ---
disp('--- Loading behavioral and FOCA data ---');
% Input data should include:
% SubjectsData: N × 400 (flattened FOCA matrices)
% MeanFD, TimeInterval, sex, scan_age, Cognitive, Language, Motor
load('./data/dHCP/BehaviorPrediction/Neonatal_FOCA_Behavior_Data.mat');

nSub = 238;     % number of neonates
nNet = 20;      % number of functional networks

SubjectsDataMatrix = reshape(SubjectsData, nSub, nNet, nNet);

% Ensure covariates are column vectors
MeanFD       = MeanFD(:);
TimeInterval = TimeInterval(:);
sex          = sex(:);
scan_age     = scan_age(:);
Cognitive    = Cognitive(:);
Language     = Language(:);
Motor        = Motor(:);

disp(['Data loaded: ' num2str(nSub) ' subjects, ' num2str(nNet) ' networks.']);

%% --- Step 2. Initialize result matrices ---
Cognitive_r_p = nan(nNet,2);
Language_r_p  = nan(nNet,2);
Motor_r_p     = nan(nNet,2);

%% --- Step 3. Regression for each network ---
disp('--- Running multiple linear regression for each network ---');
for i = 1:nNet
    % Extract FOCA profile for network i
    NetData = squeeze(SubjectsDataMatrix(:,:,i));
    NetData(NetData == 0) = 1;  % avoid zero vectors
    X = [NetData, MeanFD, TimeInterval, sex, scan_age];
    X = double(X);

    % --- Cognitive prediction ---
    stats = regstats(Cognitive, X, 'linear');
    beta = stats.beta;
    Pred_Cog = [ones(size(X,1),1), NetData] * beta(1:(nNet+1));
    [Cognitive_r_p(i,1), Cognitive_r_p(i,2)] = corr(Cognitive, Pred_Cog, 'rows','complete');

    % --- Language prediction ---
    stats = regstats(Language, X, 'linear');
    beta = stats.beta;
    Pred_Lang = [ones(size(X,1),1), NetData] * beta(1:(nNet+1));
    [Language_r_p(i,1), Language_r_p(i,2)] = corr(Language, Pred_Lang, 'rows','complete');

    % --- Motor prediction ---
    stats = regstats(Motor, X, 'linear');
    beta = stats.beta;
    Pred_Mot = [ones(size(X,1),1), NetData] * beta(1:(nNet+1));
    [Motor_r_p(i,1), Motor_r_p(i,2)] = corr(Motor, Pred_Mot, 'rows','complete');
end

%% --- Step 4. Save results ---
disp('--- Saving predictive correlation results ---');
outputDir = './results/BehaviorPrediction_GLM';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

save(fullfile(outputDir, 'Cognitive_r_p.mat'), 'Cognitive_r_p');
save(fullfile(outputDir, 'Language_r_p.mat'),  'Language_r_p');
save(fullfile(outputDir, 'Motor_r_p.mat'),     'Motor_r_p');

disp('--- Behavioral prediction analysis complete ---');
