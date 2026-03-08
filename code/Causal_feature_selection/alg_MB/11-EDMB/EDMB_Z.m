function [MB,test,time] = EDMB_Z(Data,target,alpha,ns,p,maxK)

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
%
%

if (nargin == 3)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end

start=tic;
test = 0;

sp=[];
MB=[];

spouse = cell(1, p);

% step1:find the candidate PC sets
% Use PCMB and MBOR algorithm to obtain PC sets

[PCMB_CPC,test1,~,dSep] = PCMB_PC_Z(Data,target,alpha,ns,p,maxK);
[MBOR_CPC,test2,~] = MBOR_PC_Z(Data,target,alpha,ns,p,maxK);


test = test1 + test2;

% step2:obtain the intersection,union and difference of the PC sets
PCinter = myintersect(PCMB_CPC,MBOR_CPC);
PCunion = myunion(PCMB_CPC,MBOR_CPC);
PCdiff = mysetdiff(PCunion,PCinter);

% step3:remove wrong PC nodes from PCdiff

%NonTPC euqals to U\PCunion\{T}
NonTPC = mysetdiff(1:p,myunion(PCunion,target)); 
remove = [];

for i = 1:length(PCunion)
    Y = PCunion(i);
    break_flag=0;
    for j=1:length(NonTPC)
        X = NonTPC(j);
        % S is the conditoned set or separated set
        S =  myunion(dSep{1,X}, Y);      
        test = test +1;
        % Judge independence or conditional independence among variables
        [pval]=my_fisherz_test(target,X,S,Data,ns,alpha); 
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
                   [pval]=my_fisherz_test(Y,target,S,Data,ns,alpha);      
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
PCdiff = mysetdiff(PCdiff,remove);
% Obtain PCselect
PCselect = myunion(PCinter,PCdiff);

% step4:find the spouse nodes

% some predefine matrix to store data
already_calculated_PCD=ones(1,p);
all_PCD=cell(1,p);

for i=1:length(PCselect)

    [pc_tmp,ntest2]=GetPC_Z(Data,PCselect(i),alpha, already_calculated_PCD, all_PCD, ns, p, maxK);
    %[pc_tmp,ntest2]=GetPCD_G2(Data,pc(i),alpha, ns, p, maxK);
     test=test+ntest2;
     for j=1:length(pc_tmp)

         if isempty(find(PCselect==pc_tmp(j), 1))&& pc_tmp(j)~=target && isempty(find(dSep{pc_tmp(j)}==PCselect(i), 1))

             [pval]=my_fisherz_test(pc_tmp(j),target,[dSep{pc_tmp(j)},PCselect(i)],Data,ns,alpha);

             if isnan(pval)
                 CI=0;
             else
                 if pval<=alpha
                     CI=0;
                 else
                     CI=1;
                 end
             end

             test=test+1;
             if CI==0
                 sp=myunion(sp,pc_tmp(j));
                 spouse{pc_tmp(j)} = myunion(spouse{pc_tmp(j)}, PCselect(i));
             end
         end
     end
end

% step6: 输出配偶信息到文件
fid = fopen(['EDMB_RNA_spouse_links_target' num2str(target) '.txt'], 'w');
fprintf(fid, 'Y\tA\tTarget\n');
for A = 1:p
    if ~isempty(spouse{1,A})
        for y = spouse{1,A}
            fprintf(fid, '%d\t%d\t%d\n', y, A, target);
        end
    end
end
fclose(fid);


%step5:obtain the MB set

MB=myunion(PCselect,sp);

time=toc(start);

end

