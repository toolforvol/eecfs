% 
% function [MB,test,time] = HITONMB_Z(Data,target,alpha,samples,p,maxK)
% %
% % HITONMB_Z finds the Markov blanket of target node on continuous data
% %
% % INPUT :
% %       Data is the data matrix
% %       target is the index of target node
% %       alpha is the significance level
% %       samples is the number of data samples
% %       p is the number of nodes
% %       maxK is the maximum size of conditioning set
% %
% % OUTPUT:
% %       MB is the Markov blanket of the target
% %       test is the number of conditional independence tests
% %       time is the runtime of the algorithm
% %
% %
% 
% 
% if (nargin == 3)
%    [samples,p]=size(Data);
%    maxK=3;
% end
% 
% start=tic;
% 
% test=0;
% 
% sp=[];
% MB=[];
% 
% 
% [pc,ntest1,~,sepset]=HITONPC_Z(Data,target,alpha, samples, p, maxK);
%  test=test+ntest1;
%  MB=[MB pc];
% 
% 
% for i=1:length(pc)
% 
%     [pc_tmp,ntest2]=HITONPC_Z(Data,pc(i),alpha, samples, p, maxK);
%      test=test+ntest2;
%      for j=1:length(pc_tmp)
% 
%          if isempty(find(pc==pc_tmp(j), 1))&& pc_tmp(j)~=target && isempty(find(sepset{pc_tmp(j)}==pc(i), 1))
% 
%              [CI]=my_fisherz_test(pc_tmp(j),target,[sepset{pc_tmp(j)},pc(i)],Data,samples,alpha);
% 
%              if isnan(CI)
%                  CI=0;
%              end
% 
%              test=test+1;
%              if CI==0
%                  sp=myunion(sp,pc_tmp(j));
%              end
%          end
%      end
% end
% disp(MB)
% MB=myunion(MB,sp);
% time=toc(start);

function [MB,test,time] = HITONMB_Z(Data,target,alpha,samples,p,maxK)
%
% Modified HITONMB_Z — with spouse saving logic consistent with MMMB_Z
%

if (nargin == 3)
   [samples,p] = size(Data);
   maxK = 3;
end

start = tic;

test = 0;
sp = [];          % spouse set
MB = [];          % Markov blanket

% 可选：打开日志（如不需要可注释）
logFile = fopen(['HITONMB_PROTEIN_causal_structure_log_2362_' num2str(target) '.txt'], 'w');
fprintf(logFile, "=== HITONMB Spouse Log ===\n");
fprintf(logFile, "Target Node: %d\n\n", target);

%----------------------------------------------------
% Step 1: obtain PC(target)
%----------------------------------------------------
[pc,ntest1,~,sepset] = HITONPC_Z(Data,target,alpha,samples,p,maxK);
test = test + ntest1;
MB = myunion(MB, pc);

fprintf(logFile, "[HITONPC] PC(target %d): %s\n\n", target, mat2str(pc));

%----------------------------------------------------
% Step 2: find spouses (same as MMMB_Z logic)
%----------------------------------------------------
for i = 1:length(pc)

    A = pc(i);   % 中间节点（与 MMMB_Z 完全一致）
    [pc_tmp,ntest2] = HITONPC_Z(Data, A, alpha, samples, p, maxK);
    test = test + ntest2;

    fprintf(logFile, "[HITONPC] PC(%d): %s\n", A, mat2str(pc_tmp));

    for j = 1:length(pc_tmp)

        Y = pc_tmp(j);

        % 条件1: Y 不是 target 的 PC
        if ~isempty(find(pc == Y, 1))
            continue
        end

        % 条件2: Y ≠ target
        if Y == target
            continue
        end

        % 条件3: sepset(Y) 中不能包含 A（与 MMMB_Z 完全一致）
        if ~isempty(find(sepset{Y} == A, 1))
            continue
        end

        %------------------------------------------------------------
        %   Fisher_Z test 检查 Y ⫫ target | { sepset(Y), A }
        %   逻辑与 MMMB_Z 完全相同
        %------------------------------------------------------------
        condset = [sepset{Y}, A];
        [CI] = my_fisherz_test(Y, target, condset, Data, samples, alpha);

        if isnan(CI)
            CI = 0;
        end

        test = test + 1;

        if CI == 0      % 依赖 → 认为 Y 是 target 关于 A 的配偶（与 MMMB_Z 一致）

            sp = myunion(sp, Y);

            fprintf(logFile, ...
                "  [Spouse] Node %d is a spouse of target %d via middle node %d\n", ...
                Y, target, A);
        end
    end

    fprintf(logFile, "\n");
end

%----------------------------------------------------
% Step 3: MB = PC ∪ Spouse
%----------------------------------------------------
MB = myunion(MB, sp);

time = toc(start);
fclose(logFile);

end
