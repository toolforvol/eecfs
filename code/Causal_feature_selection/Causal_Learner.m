function [Result1, Result2, Result3] = ...
    Causal_Learner(input_alg_name, data, alpha, data_type, target)

if nargin < 4
    error('Input parameters are not valid.');
end

addpath(genpath(pwd));

% Maximum size of conditioning set
maxK = 3;

% Data size
[samples, p] = size(data);

% Node size
ns = max(data);

% ************ Valid Algorithm ************
valid_algorithms = {
    'MMMB'
    'HITONMB'
    'PCMB'
    'IPCMB'
    'MBOR'
    'STMB'
    'BAMB'
    'EEMB'
    'CFS_MI'
    'EDMB'
    'CFS_MMPC'
    'EECFS'
};
if ~ismember(input_alg_name, valid_algorithms)
    error('%s is not a valid algorithm name', ...
        input_alg_name);
end
fprintf('\nMarkov blanket learning by %s\n\n', ...
    input_alg_name);

% ************ Construct Function and Run************
if strcmp(data_type,'dis')
    algorithm=str2func(strcat(input_alg_name,'_G2')); % use G2 test
    [MB, test, time] = algorithm(data,target,alpha,ns,p,maxK);
elseif strcmp(data_type,'con')
    algorithm=str2func(strcat(input_alg_name,'_Z')); % use Z test
    [MB, test, time] = algorithm (data,target,alpha,samples,p,maxK);
end

Result1 = MB;
Result2 = test;
Result3 = time;
end