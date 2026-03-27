% ************ Dataset ************ 
% | Index |   Name   |
% |  ---  |    ---   |
% |   1   |  Cancer  |
% |   2   |Earthquake|
% |   3   |  Survey  |
% |   4   |   Asia   |
% |   5   |  Sachs   |
% |   6   | Insurance|
% |   7   |   Water  |
% |   8   |  Mildew  |
% |   9   |  Barley  |
% |   10  |Hailfinder|
% |   11  |  Hepar2  |
% |   12  | Win95pts |
% |   13  |Pathfinder|
% |   14  |   Andes  |
% |   15  |   Link   |
% |   16  |   Munin  |
% data_samples=500, 1000, 1500, ..., 5000


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


clear all
clc
close all

% Index of the dataset
data_index = '1';
% Name of data
data_name = 'cancer';
% Samples of data
data_samples=5000;
% Name of algorithm
alg_name = 'MMMB';
% Significance level
alpha = 0.01;

% Index of target node. If it is global structure learning, this parameter is not needed
target = 1;   

% Path of the data set
data_path=strcat('../../data/Benchmark_NB_dataset/Result/', data_index, '-', data_name, '_', num2str(data_samples),'.txt');
if exist(data_path,'file')==0
     fprintf('\n%s does not exist.\n\n', strcat(data_path));
     return;
end

% Load data according to the path
% data needs to start from 0
data = importdata(data_path) + 1;

% Causal_Learner
[Result1, Result2, Result3]=Causal_Learner(alg_name, data, alpha, 'dis', target);
% Result1 is learned target's Markov blanket.
% Result2 is the number of conditional independence tests
% Result3 is running time


% ************ Evaluation ************
% Path of the graph
graph_path=strcat('../../data/Benchmark_NB_dataset/Result/',  data_index, '-', data_name, '_graph.txt');
if exist(graph_path, 'file')==0
     fprintf('\n%s does not exist.\n\n', graph_path);
     return;
end

% Load graph (true DAG) according to the path
graph = importdata(graph_path);

% Evaluate Markov blanket
MB=Result1;
[adj_F1,adj_precision,adj_recall]=evaluation_MB(MB, target, graph);
fprintf('\nThe learned Markov blanket of target %.0f is [', target);
for i=1:length(MB)
    if i==length(MB)
        fprintf('%d', MB(i));
    else
        fprintf('%d\t', MB(i));
    end
end
fprintf(']\n\nadj_F1=%.2f, adj_precision=%.2f, adj_recall=%.2f\n', adj_F1, adj_precision, adj_recall);
fprintf('\nThe number of conditional independence tests is %.0f.\n', Result2);
fprintf('\nElapsed time is %.2f seconds.\n\n\n', Result3);

% Save results
result_str = strcat('../../result/Benchmark_NB_dataset/', alg_name, '_metrics_', data_name, '_', data_samples, '.txt');
fid = fopen(result_str, 'wt');
fprintf(fid, 'F1\tPrecision\tRecall\tTest\tTime\n');
fprintf(fid, '%.4f\t%.4f\t%.4f\t%d\t%.4f', ...
        adj_F1, adj_precision, adj_recall, Result2, Result3);
fclose(fid);