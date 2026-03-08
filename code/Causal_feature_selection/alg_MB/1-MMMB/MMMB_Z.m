
function [MB,test,time] = MMMB_Z(Data,target,alpha,samples,p,maxK)
%
% MMMB_Z finds the Markov blanket of target node on continuous data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       samples is the number of data samples
%       p is the number of nodes
%       maxK is the maximum size of conditioning set
%
% OUTPUT:
%       MB is the Markov blanket of the target
%       test is the number of conditional independence tests
%       time is the runtime of the algorithm
%
%


if (nargin == 3)
   [samples,p]=size(Data);
   maxK=3;
end

start=tic;

test=0;

sp=[];
MB=[];

%logFile = fopen('MMMB_PROTEIN_causal_structure_log_2362.txt', 'w');  % 打开/覆盖文件
%fprintf(logFile, '=== Causal Feature Selection Log ===\n');
%fprintf(logFile, 'Target Node: %d\n', target);

[pc,ntest1,~,sepset]=MMPC_Z(Data,target,alpha, samples, p, maxK);
test=test+ntest1;
MB=[MB pc];

%fprintf(logFile, '\n[MMPC] PC of target %d: %s\n', target, mat2str(pc));
 
for i=1:length(pc)

    [pc_tmp,ntest2]=MMPC_Z(Data,pc(i),alpha, samples, p, maxK);
     test=test+ntest2;

    %fprintf(logFile, '[MMPC] PC of node %d: %s\n', pc(i), mat2str(pc_tmp));

     for j=1:length(pc_tmp)
         
         if isempty(find(pc==pc_tmp(j), 1))&& pc_tmp(j)~=target && isempty(find(sepset{pc_tmp(j)}==pc(i), 1))
             
             [CI]=my_fisherz_test(pc_tmp(j),target,[sepset{pc_tmp(j)},pc(i)],Data,samples,alpha);
             
             if isnan(CI)
                 CI=0;
             end
             
             test=test+1;
             if CI==0
                 sp=myunion(sp,pc_tmp(j));
                 %fprintf(logFile, '[Spouse] Node %d is a spouse of target %d via middle node %d\n', ...
                 %   pc_tmp(j), target, pc(i));
             end
         end
     end
end
disp(MB)

MB=myunion(MB,sp);
time=toc(start);

%fclose(logFile);  % 关闭日志文件