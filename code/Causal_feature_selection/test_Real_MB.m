% ************ Dataset ************ 
% | Index |   Name  |
% |  ---  |    ---  |
% |   1   |   DNA   |
% |   2   |   RNA   |
% |   3   | Protein |

% ************ Algorithm ************
% MMMB <
% HITONMB <
% PCMB <
% IPCMB
% MBOR
% STMB
% BAMB <
% EEMB <
% CFS_MI
% CFS_MMPC <
% EDMB <
% EECFS <

clear all
clc
close all
 
% Name of data
data_name = 'DNA';
% Samples of data
data_samples = '2362';
% Name of algorithm
alg_name = 'EECFS';
% Significance level
alpha = 0.01;
% Index of target node. If it is global structure learning, this parameter is not needed
target = 1;   

% ************ Load data ************
% data_path = strcat('../../data/', data_name, '/', data_name, '_', data_samples, '.txt');
data_path = '../../data/DNA/DNA_2362_toy.txt';
if exist(data_path, 'file')==0
     fprintf('\n%s does not exist.\n\n', data_path);
     return;
end
% data needs to start from 0
data = importdata(data_path) + 1;


% ************ Causal feature selection ************
[Result1, Result2, Result3]=Causal_Learner(alg_name, data, alpha, 'con', target);
% Markov blanket learning 
% Result1 is learned target's Markov blanket.
% Result2 is the number of conditional independence tests
% Result3 is running time

% Print the MB
MB = Result1;
fprintf('\nThe learned Markov blanket of target %.0f is [',target);
for i = 1 : length(MB)
    if i == length(MB)
        fprintf('%d', MB(i));
    else
        fprintf('%d\t', MB(i));
    end
end
fprintf(']');

% Save the MB indice
filepath = fullfile('../../result/DNA/', sprintf('%s_MB_indice.txt', alg_name));
fid = fopen(filepath, 'w');
fprintf(fid,'%d ',MB);
fclose(fid);

% ************ CITs and Runtime ************
fprintf('\nThe number of conditional independence tests is %.0f.\n',Result2);
fprintf('\nElapsed time is %.2f seconds.\n\n\n',Result3);