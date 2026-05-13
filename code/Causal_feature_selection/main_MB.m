% ************ Dataset ************ 
% | Index |     Name     |
% |  ---  |      ---     |
% |   1   |      Wine    |
% |   2   |      Heart   |
% |   3   |    Congress  |
% |   4   |      Spect   |
% |   5   |      Wdbc    |
% |   6   |     Krvskp   |
% |   7   |      Sonar   |
% |   8   |     Splice   |
% |   9   |   Bankrupty  |
% |   10  |      Seme    |
% |   11  |    Madelon   |
% |   12  |       Hiva   |
% |   13  |Ovarian-cancer|
% |   14  |   Leukemia   |
% |   15  |  Lung-cancer |
% |   16  | Breast-cancer|
% |   17  |     Dexter   |

% ************ Algorithm ************
% MMMB
% HITONMB
% PCMB
% IPCMB
% MBOR
% STMB
% BAMB
% EEMB
% CFS_MI
% CFS_MMPC
% EDMB
% EECFS

% currentFolder = fileparts(mfilename('fullpath'));
% addpath(genpath(fullfile(currentFolder, './alg_MB')));
% addpath(genpath(fullfile(currentFolder, './evaluation/eva_real_MB')));
% addpath(genpath(fullfile(currentFolder, './common')));
% addpath(genpath(fullfile(currentFolder, 'P')));
% 
% clc;
% clear;
% 
% alpha = 0.01;
% maxK = 3;
% % TODO 修改代码：
% data_name = 'Wine'; % 遍历所有数据集
% dataset_index = '1'; % 遍历所有索引
% alg_name = 'BAMB'; % 遍历所有算法
% fprintf('%s\n\n', data_name);
% 
% % ************ Load data ************ 
% filedir='../../data/Real_world_dataset/';
% data_label=load(strcat(strcat(filedir, dataset_index, '-', data_name),'/data_labels.mat'));
% cv_indices = load(strcat(strcat(filedir, dataset_index, '-', data_name),'/cv10_indices.mat'));
% 
% data=data_label.data;
% label=data_label.labels;
% indice=cv_indices.indices;
% 
% % ************ Check label balance ************ 
% size_0=length(find(label==0))/length(label);
% size_1=length(find(label==1))/length(label);
% fprintf('%.2f/%.2f\n\n',size_0,size_1);
% 
% % ************ Run causal feature selection ************ 
% disp(alg_name)
% 
% switch alg_name
%     case 'MMMB'
%         cv_nbknnsvm_MMMB(data, label, indice, data_name, alpha, maxK);
%     case 'HITONMB'
%         cv_nbknnsvm_HITONMB(data,label,indice,data_name,alpha,maxK);
%     case 'PCMB'
%         cv_nbknnsvm_PCMB(data,label,indice,data_name,alpha,maxK);
%     case 'IPCMB'
%         cv_nbknnsvm_IPCMB(data,label,indice,data_name,alpha,maxK);
%     case 'MBOR'
%         cv_nbknnsvm_MBOR(data,label,indice,data_name,alpha,maxK);
%     case 'STMB'
%         cv_nbknnsvm_STMB(data,label,indice,data_name,alpha,maxK);
%     case 'BAMB'
%         cv_nbknnsvm_BAMB(data,label,indice,data_name,alpha,maxK);
%     case 'EEMB'
%         cv_nbknnsvm_EEMB(data,label,indice,data_name,alpha,maxK);
%     case 'CFS_MI'
%         cv_nbknnsvm_CFS_MI(data,label,indice,data_name,alpha,maxK);
%     case 'CFS_MMPC'
%         cv_nbknnsvm_CFS_MMPC(data,label,indice,data_name,alpha,maxK);
%     case 'EDMB'
%         cv_nbknnsvm_EDMB(data,label,indice,data_name,alpha,maxK);
%     case 'EECFS'
%         cv_nbknnsvm_EECFS(data,label,indice,data_name,alpha,maxK);
%     otherwise
%         error('Unknown algorithm: %s', alg_name);
% end


currentFolder = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(currentFolder, './alg_MB')));
addpath(genpath(fullfile(currentFolder, './evaluation/eva_real_MB')));
addpath(genpath(fullfile(currentFolder, './common')));
addpath(genpath(fullfile(currentFolder, 'P')));

clc;
clear;

alpha = 0.01;
maxK = 3;

% =========================
% Dataset list
% =========================

% dataset_names = { ...
% 'ovarian',...
% 'Leukemia',...
% 'Lungcancer',...
% 'Breast_cancer',...
% 'Dexter'};
dataset_names = { ...
'Heart'};

dataset_indices = { ...
'2'};

% =========================
% Algorithm list
% =========================
% algorithms = { ...
% 'MMMB',...
% 'HITONMB',...
% 'PCMB',...
% 'IPCMB',...
% 'STMB',...
% 'EEMB',...
% 'CFS_MI',...
% 'CFS_MMPC',...
% 'EECFS'};
algorithms = { ...
'EECFS'};

% =========================
% CI output file (global)
% =========================

ci_global_file = '../../result/revise/All_CI_results.csv';

fid_global = fopen(ci_global_file,'wt');

fprintf(fid_global,...
'Algorithm,Dataset,Classifier,Mean,CI,Lower,Upper\n');

fclose(fid_global);

% =========================
% Main loop
% =========================

for d = 1:length(dataset_names)

    data_name = dataset_names{d};
    dataset_index = dataset_indices{d};

    fprintf('\n====================\n');
    fprintf('Dataset: %s\n', data_name);
    fprintf('====================\n');

    % ===== Load data =====

    filedir='../../data/Real_world_dataset/';

    data_label = load( ...
        strcat(filedir,...
        dataset_index,'-',data_name,...
        '/data_labels.mat'));

    cv_indices = load( ...
        strcat(filedir,...
        dataset_index,'-',data_name,...
        '/cv10_indices.mat'));

    data   = data_label.data;
    label  = data_label.labels;
    indice = cv_indices.indices;

    % ===== Label balance =====

    size_0 = length(find(label==0))/length(label);
    size_1 = length(find(label==1))/length(label);

    fprintf('Label ratio: %.2f / %.2f\n', ...
        size_0, size_1);

    % =========================
    % Loop algorithms
    % =========================

    for a = 1:length(algorithms)

        alg_name = algorithms{a};

        fprintf('\nRunning Algorithm: %s\n', ...
            alg_name);
        cv_common_nested( ...
                    data,label,indice,...
                    data_name,alpha,maxK);

        % switch alg_name
        % 
        %     case 'MMMB'
        %         cv_nbknnsvm_MMMB( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'HITONMB'
        %         cv_nbknnsvm_HITONMB( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'PCMB'
        %         cv_nbknnsvm_PCMB( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'IPCMB'
        %         cv_nbknnsvm_IPCMB( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'MBOR'
        %         cv_nbknnsvm_MBOR( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'STMB'
        %         cv_nbknnsvm_STMB( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'BAMB'
        %         cv_nbknnsvm_BAMB( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'EEMB'
        %         cv_nbknnsvm_EEMB( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'CFS_MI'
        %         cv_nbknnsvm_CFS_MI( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'CFS_MMPC'
        %         cv_nbknnsvm_CFS_MMPC( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'EDMB'
        %         cv_nbknnsvm_EDMB( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);
        % 
        %     case 'EECFS'
        %         cv_nbknnsvm_EECFS( ...
        %             data,label,indice,...
        %             data_name,alpha,maxK);

        % end

    end

end

disp('All datasets and algorithms completed.');