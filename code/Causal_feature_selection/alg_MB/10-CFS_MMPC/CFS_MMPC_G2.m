
function [MB,test,time] = CFS_MMPC_G2(Data,target,alpha,ns,p,maxK)
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



[pc,test,~,sepset] = MMPC_G2_1224(Data,target,alpha,ns,p,maxK,mysetdiff(1:p,target));
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
        [pval,dep_sp(A,Y)]=my_g2_test(Y,target,myunion(sepset{Y},A),Data,ns,alpha);    
        
           

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
            [pval]=my_g2_test(Y,A,[],Data,ns,alpha); %前面是有条件集的，现在条件集为空集判断Y和A是否之间相连，相连则认为Y是目标变量关于A的配偶
            if pval<=alpha
                spouse{1,A} = myunion(spouse{1,A},Y); %
            end

            
        end   
        
    end
    
    

    
    
end


for i =1:length(pc)
    A = pc(i);

    [spouse{1,A},test1] = MMPC_G2_1224(Data,A,alpha,ns,p,maxK,spouse{1,A}); %返回变量A的PC集合赋值给spouse{1,A}
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
            [pval]=my_g2_test(Y,target,TestSet,Data,ns,alpha);      
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


