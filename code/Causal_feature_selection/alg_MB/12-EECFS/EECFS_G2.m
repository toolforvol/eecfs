function [MB,test,time] = EECFS_G2(Data,target,alpha,ns,p,maxK)

%
% EDMB_G2 finds the Markov blanket of target node on discrete data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       ns is the size array of each node
%       p is the number of nodes
%       maxK is the maximum size of conditioning set
%
% OUTPUT:
%       MB is the Markov blanket of the target
%       test is the number of conditional independence tests
%       time is the runtime of the algorithm


if (nargin == 3) % 如果只输入前三个参数，默认ns=每个变量最大取值、p=变量数、maxK=3
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end

start=tic;
dep_sp=zeros(p,p); % 记录依赖程度
spouse = cell(1,p); % 每个变量对应的spouse集合

% step1: find the candidate PC sets using PCMB and MBOR algorithm
[PCMB_CPC, test1, ~, dSep] = PCMB_PC_G2(Data,target,alpha,ns,p,maxK); 
[MBOR_CPC, test2, ~, ~] = MBOR_PC_G2(Data,target,alpha,ns,p,maxK); %若要使用MBOR的dSep将PCSuperSet_G2.m中19行的{}改成[]即可
test = test1 + test2; % 综合 条件独立检验 （CITs） 统计

% step2: obtain the intersection, union and difference of the PC sets
PCinter = myintersect(PCMB_CPC,MBOR_CPC);
PCunion = myunion(PCMB_CPC,MBOR_CPC);
PCdiff = mysetdiff(PCunion,PCinter);

% step3: remove wrong PC nodes from PCdiff
%NonTPC euqals to U \ PCunion \ {T}, 即在U中，但不在PCunion，也不是Target的变量。
NonTPC = mysetdiff(1:p, myunion(PCunion, target));  % 1:p 代表所有变量索引, myunion表示target本身 + PC候选集合
remove = [];
for i = 1:length(PCunion) % 针对PCunion中每个元素
    Y = PCunion(i);
    break_flag=0;
    for j = 1:length(NonTPC) % 针对所有非target PC元素
        X = NonTPC(j);
        % S is the conditoned set or separated set S是条件集或分离集
        S = myunion(dSep{1,X}, Y);   
        test = test + 1;
        % Judge independence or conditional independence among variables
        [pval]=my_g2_test(target, X, S, Data, ns, alpha); 
        if isnan(pval)        
            CI=0;
        else
            if pval<=alpha
                CI=0;
            else
                CI=1;
            end
        end
        % According to independence result,remove wrong PC nodes from PCdiff
        if(CI==0)    
           SZ = mysetdiff( myunion(PCunion,X),Y );
           cutSetSize = 0;
           while length(SZ) >= cutSetSize&&cutSetSize<=maxK
               SS = subsets1(SZ, cutSetSize);   
               for si=1:length(SS)
                   S = SS{si};
                   test=test+1;
                   [pval]=my_g2_test(Y,target,S,Data,ns,alpha); 
                   if isnan(pval)      
                       CI=0;
                   else
                       if pval<=alpha
                           CI=0;
                       else
                           CI=1;
                       end
                   end
                   if(CI==1)  
                       remove = myunion(remove,Y);
                       break_flag=1; 
                       break;
                   end
               end
               if( break_flag==1 )
                   break;
               end
               cutSetSize = cutSetSize + 1;       
           end
           if( break_flag==1 )
               break;
           end
        end
    end
end
% Process PCdiff set
PCdiff = mysetdiff(PCdiff, remove);
% Obtain PCselect
PCselect = myunion(PCinter,PCdiff);

% step4:find the spouse nodes(different from EDMB)
NoPC = mysetdiff(1:p,[PCselect,target]);

for i =1:length(PCselect)
    % A is the current node in the target variable PC set
    A = PCselect(i); % A equals to X in the CFS framework of the literature
    for j=1:length(NoPC)
        Y = NoPC(j); % Y is the current node in the non PC set of the target

        if ~isempty(find(dSep{Y}==A, 1))%之前为dSep{Y}==A
            continue
        end 

        test = test+1;
        [pval,dep_sp(A,Y)]=my_g2_test(Y,target,myunion(dSep{Y},A),Data,ns,alpha);   % before dSep{Y}
        if isnan(pval)
            CI=1;
        else
            if pval<=alpha
                CI=0;
            else
                CI=1;
            end
        end

        if CI==0 % If it is a dependency, then Y is a potential spouse of T with respect to X, and X will become a variable with multiple parent nodes
            test = test+1;
            [pval]=my_g2_test(Y,A,[],Data,ns,alpha); % Previously, there was a conditional set, but now the condition set is an empty set to determine whether Y and A are connected. If they are connected, it is considered that Y is the spouse of the target variable with respect to A
            if pval<=alpha
                spouse{1,A} = myunion(spouse{1,A},Y);
            end
            
        end   
        
    end
    
end


for i =1:length(PCselect)
    A = PCselect(i);

    [spouse{1,A},test1] = MMPC_G2_1224(Data,A,alpha,ns,p,maxK,spouse{1,A}); % Return the PC collection of variable A and assign it to spouse{1,A}
        
    test = test+test1;
    % This step is equivalent to resetting spot {1, A}, initially only obtaining a potential set of spouses
    % As the input for the second learning of the PC set, the final spouse set that needs to be filtered is obtained at this time
    % Equivalent to updating the spouse set
 %--------------------------------------2---------------
    tmp_PC = spouse{1,A};
    for j=(length(spouse{1,A})):-1:1
        break_flag=0;
%         j
%         disp(length(spouse{1,A}));
        X = spouse{1,A}(j);

        CanPC=mysetdiff(myunion(tmp_PC,PCselect), [X,A]);

        cutSetSize = 1;

        while length(CanPC) >= cutSetSize&&cutSetSize+1<=maxK

            SS = subsets1(CanPC, cutSetSize);    

            for si=1:length(SS)
                Z = SS{si};

                test=test+1;

                [pval]=my_g2_test(target,X,myunion(Z,A),Data,ns,alpha);    
                
                if isnan(pval)
                    CI=0;
                else
                    if pval<=alpha
                        CI=0;
                    else
                        CI=1;
                    end
                end

                if (CI==1)      
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
    
    spouse{1,A}=tmp_PC;% Retain the final set of spouses
 
end

fid = fopen(['EECFS_DNAdis_spouse_links_target' num2str(target) '.txt'], 'w');
fprintf(fid, 'Y\tA\tTarget\n');
for A = 1:p
    if ~isempty(spouse{1,A})
        for y = spouse{1,A}
            fprintf(fid, '%d\t%d\t%d\n', y, A, target);
        end
    end
end
fclose(fid);



% %The final step in the framework is to remove false positive nodes
% pc_tmp=PCselect;
% 
% 
% for i=1:length(PCselect) 
%     Y = PCselect(i);
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
%             [pval]=my_g2_test(Y,target,TestSet,Data,ns,alpha);      
%             if isnan(pval)
%                 CI=0;
%             else
%                 if pval<=alpha
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
% PCselect=pc_tmp; 
     
MB=myunion(PCselect,cell2mat(spouse));


time=toc(start);

end

