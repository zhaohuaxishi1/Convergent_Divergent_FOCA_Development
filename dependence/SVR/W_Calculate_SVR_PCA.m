function [w_Brain, model_All] = W_Calculate_SVR_PCA(Subjects_Data, Subjects_Scores, Pre_Method, C_Parameter, Netnum,Netfea,Voxelfea,Indicatorfea,Reduction_Method,randNet,ResultantFolder)
%
% Subject_Data:
%           m*n matrix
%           m is the number of subjects
%           n is the number of features
%
% Subject_Scores:
%           the continuous variable to be predicted,[1*m]
%
% Pre_Method:
%          'Normalize', 'Scale', 'None'
%
% C_Parameter:
%          We generally use 1 as default C parameter.
%          See Cui et al., 2017, Human Brain Mapping
%
% ResultantFolder:
%          the path of folder storing resultant files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Zaixu Cui: zaixucui@gmail.com;
%                       Zaixu.Cui@pennmedicine.upenn.edu
%
% If you use this code, please cite:
%                       Cui et al., 2018, Cerebral Cortex;
%                       Cui and Gong, 2018, NeuroImage;
%                       Cui et al., 2016, Human Brain Mapping.
% (google scholar: https://scholar.google.com.hk/citations?user=j7amdXoAAAAJ&hl=zh-TW&oi=ao)
%

if nargin >= 8
    if ~exist(ResultantFolder, 'dir')
        mkdir(ResultantFolder);
    end
end

% data reduction
Subjects_Data_pca = [];
coeff_pca =  [];
for m = 1:Netnum
    % Traing test
    Subjects_Data_voxel = Subjects_Data(:,(randNet(m)-1)*Netfea+1:randNet(m)*Netfea);
    if strcmp(Reduction_Method, 'Pca')
        [coeff,score,latent,tsquared,explained,~] = pca(Subjects_Data_voxel);
        Subjects_Data_voxel_pca = Subjects_Data_voxel * coeff(:,:);
        Subjects_Data_pca = [Subjects_Data_pca,Subjects_Data_voxel_pca];
        coeff_pca =  [coeff_pca;coeff];
    end
end

NetfeaPCA = size(Subjects_Data,1) + Indicatorfea - 1;
[~, Features_Quantity] = size(Subjects_Data_voxel);


if strcmp(Pre_Method, 'Normalize')
    %Normalizing
    MeanValue = mean(Subjects_Data_pca);
    StandardDeviation = std(Subjects_Data_pca);
    [~, columns_quantity] = size(Subjects_Data_pca);
    for j = 1:columns_quantity
        Subjects_Data_pca(:, j) = (Subjects_Data_pca(:, j) - MeanValue(j)) / StandardDeviation(j);
    end
elseif strcmp(Pre_Method, 'Scale')
    % Scaling to [0 1]
    MinValue = min(Subjects_Data_pca);
    MaxValue = max(Subjects_Data_pca);
    [~, columns_quantity] = size(Subjects_Data_pca);
    for j = 1:columns_quantity
        Subjects_Data_pca(:, j) = (Subjects_Data_pca(:, j) - MinValue(j)) / (MaxValue(j) - MinValue(j));
    end
end
% bracause Subjects data contain zero value of some columns, thus will have
% NAN value after divide zero
Subjects_Data_pca(isnan(Subjects_Data_pca)==1) = 0;
% SVR
Subjects_Data_pca = double(Subjects_Data_pca);

model_All = svmtrain(Subjects_Scores, Subjects_Data_pca,['-s 3 -t 0 -c ' num2str(C_Parameter)]);


% Select coeff
coeff_pca_net = [];
for k = 1:Netnum
    coeff_pca_net = [coeff_pca_net;coeff_pca((k-1)*Voxelfea+1:k*Voxelfea,:)];
end

w_Brain_PCA = zeros(1, size(Subjects_Data_pca,2));
for m = 1 : model_All.totalSV
    w_Brain_PCA = w_Brain_PCA + model_All.sv_coef(m) * model_All.SVs(m, :);
end

w_Brain = [];
for m =1:Netnum
    w_Brain_tmp = w_Brain_PCA(1,(m-1)*NetfeaPCA+1:m*NetfeaPCA) * coeff_pca_net((m-1)*Voxelfea+1:m*Voxelfea,:)';
    w_Brain = [w_Brain,w_Brain_tmp];
end
w_Brain = w_Brain / norm(w_Brain);

if nargin >= 5
    save([ResultantFolder filesep 'w_Brain.mat'], 'w_Brain');
    save([ResultantFolder filesep 'w_Brain_PCA.mat'], 'w_Brain_PCA');
end
