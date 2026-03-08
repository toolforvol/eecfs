
function [MB,test,time] = CFS_MMPC_Z(Data,target,alpha,ns,p,maxK)
%



if (nargin == 3)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end


start=tic;

cpc=[];
pc=[];

dep=zeros(1,p);

spouse = cell(1,p);

dep_sp=zeros(p,p);

% Remove the nodes independent of the target node conditioning on
% empty set, and sort the dependent nodes


sepset=cell(1,p); %create p arrays to store every feature's spouse
test=0;

U = mysetdiff(1:p,target); %除目标变量外的所有特征索引

score=1000000;

last_added=-1;

pval=ones(1,p)*score;
dep_tmp1 = ones(1,p)*score;  %X到T的依赖程度最小值

NoPC = [];



[pc,test,~,sepset] = MMPC_Z_1224(Data,target,alpha,ns,p,maxK,mysetdiff(1:p,target));
% pc
% NoPC
% U
% NoPC = myunion(NoPC,U);
% NoPC
NoPC = mysetdiff(1:p,[pc,target]);

for i =1:length(pc)
    % A也就是目标变量PC集中的当前节点
    A = pc(i); % A equals to X in the CFS framework of the literature
    for j=1:length(NoPC)
        Y = NoPC(j); %Y是target的非PC集中的当前节点
        if ~isempty(find(sepset{Y}==A, 1))
            continue
        end 
        test = test+1;
        [pval,dep_sp(A,Y)]=my_fisherz_test(Y,target,myunion(sepset{Y},A),Data,ns,alpha);    



        if isnan(pval)
            CI=1;
        else
            if pval<=alpha
                CI=0;
            else
                CI=1;
            end
        end

        if CI==0 %如果是依赖，则Y是T关于X的潜在配偶，那么X将会成为有多个父节点的变量
            test = test+1;
            [pval]=my_fisherz_test(Y,A,[],Data,ns,alpha); %前面是有条件集的，现在条件集为空集判断Y和A是否之间相连，相连则认为Y是目标变量关于A的配偶
            if pval<=alpha
                spouse{1,A} = myunion(spouse{1,A},Y); %
            end


        end   

    end





end


for i =1:length(pc)
    A = pc(i);

    [spouse{1,A},test1] = MMPC_Z_1224(Data,A,alpha,ns,p,maxK,spouse{1,A}); %返回变量A的PC集合赋值给spouse{1,A}
    test = test+test1;
    %这一步相当于把spouse{1,A}进行重置，最开始只是得到一个潜在的配偶集合，作为第二次学习PC集合的输入，此时得到的才是最终需要进行筛选的配偶集合
    %相当于将配偶集进行更新
 %--------------------------------------2---------------






    tmp_PC = spouse{1,A};
    for j=(length(spouse{1,A})):-1:1
         break_flag=0;
%         j
%         disp(length(spouse{1,A}));
        X = spouse{1,A}(j);

        CanPC=mysetdiff(myunion(tmp_PC,pc), [X,A]);

        cutSetSize = 1;

        while length(CanPC) >= cutSetSize&&cutSetSize+1<=maxK

            SS = subsets1(CanPC, cutSetSize);    

            for si=1:length(SS)
                Z = SS{si};

                test=test+1;

                [pval]=my_fisherz_test(target,X,myunion(Z,A),Data,ns,alpha);    


                if isnan(pval)
                    CI=0;
                else
                    if pval<=alpha
                        CI=0;
                    else
                        CI=1;
                    end
                end

                if CI==1          
                    tmp_PC = mysetdiff(tmp_PC,X);
%                     sepset{X}=Z;

                    break_flag=1;
                    break;
                end
            end

            if break_flag
                break;
            end
            cutSetSize = cutSetSize + 1;
        end   

    end

    spouse{1,A}=tmp_PC;% 保留最终的配偶集



end













%框架中最后一步移除假阳性节点
pc_tmp=pc;


for i=1:length(pc) 
    Y = pc(i);

    CanPC=mysetdiff(pc_tmp, Y);

    cutSetSize = 1;            
    other_PC_break_flag=0;
    while length(CanPC) >= cutSetSize && cutSetSize<=maxK
        SS = subsets1(CanPC, cutSetSize);   
        for si=1:length(SS)
            Z = SS{si};        

            spouse_test=[];
            for k=1:length(Z)
                pc_var = Z(k);
                spouse_test = myunion(spouse_test,spouse{1,pc_var});
            end
            TestSet = myunion(Z,spouse_test);

            test=test+1;
            [pval]=my_fisherz_test(Y,target,TestSet,Data,ns,alpha);      
            if isnan(pval)
                CI=0;
            else
                if pval<=alpha
                    CI=0;
                else
                    CI=1;
                end
            end

            if CI==1          
                pc_tmp=mysetdiff( pc_tmp,Y );
                spouse{1,Y} = [];
                other_PC_break_flag=1;   
                break;
            end
        end
        if other_PC_break_flag==1
            break;
        end
        cutSetSize = cutSetSize + 1;
    end

end
pc=pc_tmp; 


MB=myunion(pc,cell2mat(spouse));


time=toc(start);



% function [MB,test,time] = CFS_MMPC_Z(Data,target,alpha,ns,p,maxK)
% %
% % Modified CFS_MMPC_Z that records spouse nodes (with logging) similarly to
% % the provided MMMB_Z implementation. The function keeps the original
% % algorithmic flow but writes events to a log file and saves spouse
% % assignments per-PC node. Final MB is union of pc and all spouse nodes.
% %
% % INPUT/OUTPUT same as original CFS_MMPC_Z
% %
% if (nargin == 3)
%    ns=max(Data);
%    [~,p]=size(Data);
%    maxK=3;
% end
% 
% start=tic;
% 
% cpc=[];
% pc=[];
% 
% dep=zeros(1,p);
% 
% spouse = cell(1,p);
% 
% dep_sp=zeros(p,p);
% 
% % Remove the nodes independent of the target node conditioning on
% % empty set, and sort the dependent nodes
% 
% sepset=cell(1,p); %create p arrays to store every feature's spouse
% test=0;
% 
% U = mysetdiff(1:p,target); %除目标变量外的所有特征索引
% 
% score=1000000;
% 
% last_added=-1;
% 
% pval=ones(1,p)*score;
% dep_tmp1 = ones(1,p)*score;  %X到T的依赖程度最小值
% 
% NoPC = [];
% 
% % --- open log file (overwrites to start fresh, similar to MMMB_Z) ---
% logFile = fopen('CFS_MMPC_DNA_causal_structure_log_2362.txt', 'w');
% if logFile~=-1
%     fprintf(logFile, '=== Causal Feature Selection Log (CFS_MMPC_Z) ===\n');
%     fprintf(logFile, 'Target Node: %d\n', target);
%     fprintf(logFile, 'alpha = %g, ns = %d, p = %d, maxK = %d\n\n', alpha, ns, p, maxK);
% else
%     warning('Could not open log file for writing. Continuing without file logging.');
% end
% 
% [pc,test,~,sepset] = MMPC_Z_1224(Data,target,alpha,ns,p,maxK,mysetdiff(1:p,target));
% fprintf(logFile, '[CFS_MMPC] PC of target %d: %s\n', target, mat2str(pc));
% 
% NoPC = mysetdiff(1:p,[pc,target]);
% 
% % --- open spouse links file to immediately record discovered spouse relations ---
% spouseLinksFile = fopen('CFS_MMPC_DNA_Spouse_Links.txt','w');
% if spouseLinksFile~=-1
%     fprintf(spouseLinksFile,'Y  A   Target'); % header (optional)
% else
%     warning('Could not open spouse links file for writing.');
% end
% 
% % First pass: find candidate spouses by testing nodes in NoPC against target
% for i =1:length(pc)
%     % A也就是目标变量PC集中的当前节点
%     A = pc(i); % A equals to X in the CFS framework of the literature
%     for j=1:length(NoPC)
%         Y = NoPC(j); %Y是target的非PC集中的当前节点
%         if ~isempty(find(sepset{Y}==A, 1))
%             continue
%         end 
%         test = test+1;
%         [pval_tmp,dep_sp(A,Y)]=my_fisherz_test(Y,target,myunion(sepset{Y},A),Data,ns,alpha);    
% 
%         if isnan(pval_tmp)
%             CI=1; % treat NaN conservatively as independent
%         else
%             if pval_tmp<=alpha
%                 CI=0;
%             else
%                 CI=1;
%             end
%         end
% 
%         if CI==0 %如果是依赖，则Y是T关于X的潜在配偶，那么X将会成为有多个父节点的变量
%             test = test+1;
%             [pval2]=my_fisherz_test(Y,A,[],Data,ns,alpha); %判断Y和A是否之间相连
%             if ~isnan(pval2) && pval2<=alpha
%                 spouse{1,A} = myunion(spouse{1,A},Y); % record candidate spouse for A
%                 % Immediately write discovered spouse link to file so it won't be lost by subsequent refinement
%                 if exist('spouseLinksFile','var') && spouseLinksFile~=-1
%                     fprintf(spouseLinksFile,'%d	%d	%d',Y,A,target);
%                 end
%                 % Log the detection (similar to MMMB_Z style)
%                 if logFile~=-1
%                     fprintf(logFile, '[CandidateSpouse] Node %d is a candidate spouse of target %d via middle node %d (pval(Y,T|...)=%.5g, pval(Y,A)=%.5g)\n', Y, target, A, pval_tmp, pval2);
%                 end
%             end
% 
%         end   
% 
%     end
% end
% 
% % Second pass: refine spouse sets by computing PC of A using initial spouse candidates
% for i =1:length(pc)
%     A = pc(i);
% 
%     % Use MMPC_Z_1224 with current spouse candidates as "initial set"
%     [sp_set,test1] = MMPC_Z_1224(Data,A,alpha,ns,p,maxK,spouse{1,A}); %返回变量A的PC集合赋值给spouse{1,A}
%     test = test+test1;
% 
%     % Log the intermediate PC returned for A
%     if logFile~=-1
%         fprintf(logFile, '[CFS_MMPC] PC learned for node %d (used initial candidates): %s\n', A, mat2str(sp_set));
%     end
% 
%     tmp_PC = sp_set;
% 
%     % prune tmp_PC: for each candidate X in spouse{1,A}, test whether target indep of X given (Z + A)
%     for j=(length(sp_set)):-1:1
%         break_flag=0;
%         X = sp_set(j);
% 
%         CanPC=mysetdiff(myunion(tmp_PC,pc), [X,A]);
% 
%         cutSetSize = 1;
% 
%         while length(CanPC) >= cutSetSize && cutSetSize+1<=maxK
% 
%             SS = subsets1(CanPC, cutSetSize);    
% 
%             for si=1:length(SS)
%                 Z = SS{si};
% 
%                 test=test+1;
% 
%                 [pval3]=my_fisherz_test(target,X,myunion(Z,A),Data,ns,alpha);    
% 
%                 if isnan(pval3)
%                     CI=0; % treat NaN as dependent (to be conservative in pruning)
%                 else
%                     if pval3<=alpha
%                         CI=0;
%                     else
%                         CI=1;
%                     end
%                 end
% 
%                 if CI==1          
%                     tmp_PC = mysetdiff(tmp_PC,X);
%                     % record the separating set if desired: sepset{X}=Z;
% 
%                     break_flag=1;
%                     if logFile~=-1
%                         fprintf(logFile, '[PruneSpouse] Removed node %d from spouse candidates of node %d; separating set: %s\n', X, A, mat2str(Z));
%                     end
%                     break;
%                 end
%             end
% 
%             if break_flag
%                 break;
%             end
%             cutSetSize = cutSetSize + 1;
%         end   
% 
%     end
% 
%     spouse{1,A}=tmp_PC; % 保留最终的配偶集
% 
%     if logFile~=-1
%         fprintf(logFile, '[FinalSpouse] Final spouse set for node %d: %s\n\n', A, mat2str(spouse{1,A}));
%     end
% 
% end
% 
% % 框架中最后一步移除假阳性节点（对pc进行最终筛选）
% pc_tmp=pc;
% 
% for i=1:length(pc) 
%     Y = pc(i);
% 
%     CanPC=mysetdiff(pc_tmp, Y);
% 
%     cutSetSize = 1;            
%     other_PC_break_flag=0;
%     while length(CanPC) >= cutSetSize && cutSetSize<=maxK
%         SS = subsets1(CanPC, cutSetSize);   
%         for si=1:length(SS)
%             Z = SS{si};        
% 
%             spouse_test=[];
%             for k=1:length(Z)
%                 pc_var = Z(k);
%                 spouse_test = myunion(spouse_test,spouse{1,pc_var});
%             end
%             TestSet = myunion(Z,spouse_test);
% 
%             test=test+1;
%             [pval4]=my_fisherz_test(Y,target,TestSet,Data,ns,alpha);      
%             if isnan(pval4)
%                 CI=0;
%             else
%                 if pval4<=alpha
%                     CI=0;
%                 else
%                     CI=1;
%                 end
%             end
% 
%             if CI==1          
%                 pc_tmp=mysetdiff( pc_tmp,Y );
%                 spouse{1,Y} = [];
%                 other_PC_break_flag=1;   
%                 if logFile~=-1
%                     fprintf(logFile, '[RemovePC] Removed PC node %d from target %d; TestSet used: %s\n', Y, target, mat2str(TestSet));
%                 end
%                 break;
%             end
%         end
%         if other_PC_break_flag==1
%             break;
%         end
%         cutSetSize = cutSetSize + 1;
%     end
% 
% end
% pc=pc_tmp; 
% 
% 
% % === After final PC & spouse extraction, write final spouse links ===
% if exist('spouseLinksFile','var') && spouseLinksFile~=-1
%     % Clear previous content & rewrite final only
%     fclose(spouseLinksFile);
%     spouseLinksFile = fopen('CFS_MMPC_DNA_Spouse_Links.txt','w');
%     fprintf(spouseLinksFile,'Y\tA\tTarget\n');
% 
%     for A=pc % final PC nodes only
%         if ~isempty(spouse{A})
%             for idx=1:length(spouse{A})
%                 Y = spouse{A}(idx);
%                 fprintf(spouseLinksFile,'%d\t%d\t%d\n',Y,A,target);
%             end
%         end
%     end
%     fclose(spouseLinksFile);
% end
% 
% % MB is union of pc and all spouse members
% all_spouses = [];
% for i=1:length(spouse)
%     if ~isempty(spouse{i})
%         all_spouses = myunion(all_spouses, spouse{i});
%     end
% end
% 
% % spouse links are recorded immediately during discovery above; ensure file is closed below before finishing.
% 
% % close spouse links file if open
% % if exist('spouseLinksFile','var') && spouseLinksFile~=-1
% %     fclose(spouseLinksFile);
% % end
% 
% MB=myunion(pc,all_spouses);
% 
% if logFile~=-1
%     fprintf(logFile, '\n[Result] Final PC set for target %d: %s\n', target, mat2str(pc));
%     fprintf(logFile, '[Result] Final spouse nodes combined: %s\n', mat2str(all_spouses));
%     fprintf(logFile, '[Result] Final MB: %s\n', mat2str(MB));
%     fclose(logFile);
% end
% 
% 
% time=toc(start);
% 
% end


