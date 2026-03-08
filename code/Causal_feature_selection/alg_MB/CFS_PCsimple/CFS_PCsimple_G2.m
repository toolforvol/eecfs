function [MB,test,time] = CFS_PCsimple_G2(Data,target,alpha,ns,p,maxK)

if (nargin == 3)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end

start=tic;

spouse = cell(1,p);

dep_sp=zeros(p,p);

% Remove nodes that are independent of the target node given an empty set, 
% and sort the dependent nodes
score=1000000;

ADJT = mysetdiff(1:p,target);

% Obtain the PC set and separation set of the target variable
[pc,test,~,sepset] = PCsimple_G2_1224(Data,target,alpha,ns,p,maxK,mysetdiff(1:p,target));     
NoPC = mysetdiff(1:p,[pc,target]); % Obtain the non-PC set after obtaining the PC set
% pc

% Loop through each node in the PC set
for i =1:length(pc)
    X = pc(i); % X is the current node in the PC set of the target
    for j=1:length(NoPC)
        Y = NoPC(j); % Y is the current node in the non-PC set of the target
        if ~isempty(find(sepset{Y}==X, 1)) % If X is in the separation set of Y, skip this pair
            continue;
        end      
              
        test = test+1;
        % Test the relationship between the target variable and Y conditioned on X and their separation set
        [pval,dep_sp(X,Y)]=my_g2_test(Y,target,myunion(sepset{Y},X),Data,ns,alpha);    

        if isnan(pval)
            CI=0;
        else
            if pval<=alpha
                CI=0;
            else
                CI=1;
            end
        end

        if CI==0     % If dependent, Y is a potential spouse of target w.r.t X

            spouse{1,X} = myunion(spouse{1,X},Y); % Add Y as a spouse of the target with respect to X
        end   
        
    end
    
    Y = pc(i); % Now Y is the current node in the target's PC set (original X)

    % After the above, Y (i.e., X) may have multiple parents
    % Use the PC set of nodes with multiple parents to identify spouses
    [spouse{1,Y},test1] = PCsimple_G2_1224(Data,Y,alpha,ns,p,maxK,spouse{1,Y}); 
    test = test+test1;
    % This step updates spouse{1,X} to the final spouse set after the second PC learning

    ADJT_length = length(spouse{1,Y});
    ADJT = spouse{1,Y}; % Y's PC set
    tmp_ADJT = ADJT;
    cutSetSize = 0;
    % Maximum conditioning set size = maxK
    while length(myunion(ADJT,pc)) > cutSetSize && cutSetSize+1 <= maxK
        for j=1:length(ADJT)
            X = ADJT(j); % Current node in Y's PC set
            % Neighbors excluding X and Y from the union of ADJT and PC
            nbrs = mysetdiff(myunion(ADJT,pc), [X,Y]); 

            SS = subsets1(nbrs, cutSetSize);    % Enumerate subsets of size cutSetSize
            for si=1:length(SS)
                S = SS{si};
                test=test+1;
                [pval]=my_g2_test(target,X,myunion(S,Y),Data,ns,alpha);     
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
                    tmp_ADJT = mysetdiff(tmp_ADJT,X); % If independent, remove X from the current set
                    ADJT_length = ADJT_length-1;
                    break; % no need to test other subsets
                end
            end

        end
        ADJT=tmp_ADJT;
        spouse{1,Y} = ADJT; % Keep the final spouse set
        cutSetSize = cutSetSize + 1; % Increase conditioning set size
    end    
end

% --- Remove false positives
pc_tmp=pc;
for i=1:length(pc) 
    Y = pc(i);

    CanPC=mysetdiff(pc, Y);

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