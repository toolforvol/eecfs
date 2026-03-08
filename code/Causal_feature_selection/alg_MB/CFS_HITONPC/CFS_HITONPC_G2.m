
function [MB,test,time] = CFS_HITONPC_G2(Data,target,alpha,ns,p,maxK)
%



if (nargin == 3)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end


start=tic;


test=0;

cpc=[];
pc=[];

dep=zeros(1,p);
sepset=cell(1,p);

spouse = cell(1,p);

dep_sp=zeros(p,p);

% Remove the nodes independent of the target node conditioning on
% empty set, and sort the dependent nodes

for i=1:p
    if i~=target
        
        test=test+1;
        [pval, dep(i)]=my_g2_test(i,target,[],Data,ns,alpha);

        if isnan(pval)
            CI=1;      
        else
            if pval<=alpha
                CI=0;
            else
                CI=1;
            end
        end
        
        if CI==1
            continue;
        end
        if CI==0
            cpc=[cpc,i];
        end
    end
end

[~,var_index]=sort(dep(cpc),'descend');

CI_dependent=[];
for i=1:length(cpc)
    CI_dependent=[CI_dependent cpc(var_index(i))];
end

% NonTPC = mysetdiff(1:p,myunion(CI_dependent,target));


% ---------------------------------------------- ADDTrue









for i=1:length(cpc)
    pc=[pc cpc(var_index(i))];
    pc_length=length(pc);
    pc_tmp=pc;
    last_break_flag=0;

    % remove false PC
    
    for j=pc_length:-1:1
        Y = pc(j);
        CanPC=mysetdiff(pc_tmp, Y);
                    
        cutSetSize = 1;             
        other_break_flag=0;
        while length(CanPC) >= cutSetSize && cutSetSize<=maxK
            SS = subsets1(CanPC, cutSetSize);    
            for si=1:length(SS)
                Z = SS{si};
                
                if cpc(var_index(i))~=Y       
                    if isempty(find(Z==cpc(var_index(i)), 1))
                        continue;
                    end
                end                 

                test=test+1;
                [pval]=my_g2_test(Y,target,Z,Data,ns,alpha);     
                if isnan(pval)
                    CI=0;
                else
                    if pval<=alpha
                        CI=0;
                    else
                        CI=1;
                    end
                end
                
                % find spouses with regard to the PC without Y 
                
                if CI==1 
                    pc_tmp=CanPC;
                    sepset{Y}=Z;

%                     NonTPC=myunion(NonTPC,Y);               
                    
                    if cpc(var_index(i))==Y
                        last_break_flag=1;
                    end
                    
                    spouse{1,Y}=[];        
                    
                    other_break_flag=1;    
                    break;
                end 
            end

            if last_break_flag||other_break_flag
                break;
            end

            cutSetSize = cutSetSize + 1;
        end

        if last_break_flag
            break;
        end
        
        % the new added node is belong to PC, 
        % find spouses with regard to this node

    end
    pc=pc_tmp;
end

NonTPC = mysetdiff(1:p,[pc,target]);

% pc
% NonTPC

M=pc;

for i =1:length(cpc)
    if isempty(find(M==cpc(var_index(i)), 1))
        continue;
    end
    A = cpc(var_index(i));
    for j=1:length(NonTPC)
        Y = NonTPC(j);

        
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

        if CI==0          
            spouse{1,A} = myunion(spouse{1,A},Y);
        end   
        
        
    end
    
    Y=cpc(var_index(i));
    
   [~,SP_index]=sort(dep_sp(Y,spouse{1,Y}),'descend'); 
    SP=[];
    for f=1:length(spouse{1,Y})
        SP=[SP spouse{1,Y}(SP_index(f))];
        SP_length=length(SP);
        SP_tmp=SP;
        SP_break_flag2=0;
        for c=SP_length:-1:1
            X = SP(c);
            CanSP=mysetdiff(SP_tmp, X);
            cutSetSize = 0;
            
            other_SP_break_flag=0;
            while length(CanSP) >= cutSetSize && cutSetSize<=maxK
                SS = subsets1(CanSP, cutSetSize);   
                for si=1:length(SS)
                    Z = SS{si};
                    
                    condset=Z;
                    
                    if spouse{1,Y}(SP_index(f))~=X
                        if isempty(find(Z==spouse{1,Y}(SP_index(f)), 1))
                            continue;
                        end
                    end
                    
                    test=test+1;
%                     X
%                     Y
%                     condset
                    [pval]=my_g2_test(X,Y,condset,Data,ns,alpha);
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
                        
                        SP_tmp=mysetdiff(SP_tmp, X);
                        
                        if spouse{1,Y}(SP_index(f))==X
                            SP_break_flag2=1;
                        end
                        other_SP_break_flag=1;
                        break;
                    end
                end
                if( SP_break_flag2==1 )
                    break;
                end
                if( other_SP_break_flag==1 )
                    break;
                end
                cutSetSize = cutSetSize + 1;
            end
            if( SP_break_flag2==1 )      
                break;
            end
        end
        SP=SP_tmp;
    end
    spouse{1,Y} = SP;
    
    
%     Phase I-3
    [~,SP_index]=sort(dep_sp(Y,spouse{1,Y}),'descend');
    SP=[];
    for f=1:length(spouse{1,Y})
        SP=[SP spouse{1,Y}(SP_index(f))];
        SP_length=length(SP);
        SP_tmp=SP;
        SP_break_flag=0;
        for c=SP_length:-1:1
            X = SP(c);
            CanSP=mysetdiff(myunion(pc,SP_tmp), [X,Y]);
            cutSetSize = 1;            
            
            other_SP_break_flag=0;
            while length(CanSP) >= cutSetSize&& cutSetSize+1<=maxK
                SS = subsets1(CanSP, cutSetSize);  
                for si=1:length(SS)
                    Z = SS{si};
                    condset=unique([Z,Y]);
                    
                    
                    if (length(intersect(pc,condset))==1)
                        continue;
                    end                    
                    
                    if spouse{1,Y}(SP_index(f))~=X    
                        if isempty(find(Z==spouse{1,Y}(SP_index(f)), 1))
                            continue;
                        end
                    end
                    
                    test=test+1;
                    [pval]=my_g2_test(X,target,condset,Data,ns,alpha);   
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
                        
                        SP_tmp=mysetdiff(SP_tmp, X);
                        
                        if spouse{1,Y}(SP_index(f))==X
                            SP_break_flag=1;
                        end
                        other_SP_break_flag=1;
                        break;
                    end
                end
                if( SP_break_flag==1 )
                    break;
                end
                if( other_SP_break_flag==1 )
                    break;
                end
                cutSetSize = cutSetSize + 1;
            end
            
            if( SP_break_flag==1 )     
                break;
            end
        end
        SP=SP_tmp;
    end
    spouse{1,Y} = SP;

      
    
    
end



% ---------------------------------------------- RMFalse
% Phase I: Remove false positives from SPST




% Phase II: Remove false positives from PCST


M=pc;
pc=[];
del_spouse=[];


for i=1:length(cpc) 
    if isempty(find(M==cpc(var_index(i)), 1))
        continue;
    end
    pc=[pc cpc(var_index(i))];
    pc_length=length(pc);
    pc_tmp=pc;
    last_PC_break_flag=0;

    for j=pc_length:-1:1
        Y = pc(j);
        CanPC=mysetdiff(pc_tmp, Y);
                    
        cutSetSize = 1;            
        other_PC_break_flag=0;
        while length(CanPC) >= cutSetSize && cutSetSize<=maxK
            SS = subsets1(CanPC, cutSetSize);   
            for si=1:length(SS)
                Z = SS{si};
                
                if cpc(var_index(i))~=Y            
                    if isempty(find(Z==cpc(var_index(i)), 1))
                        continue;
                    end
                end 
                
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
                    pc=mysetdiff( pc,Y );    %alter
                    del_spouse = spouse{1,Y};
                    spouse{1,Y} = [];
                    if cpc(var_index(i))==Y
                        last_PC_break_flag=1;   
                    end
                    other_PC_break_flag=1;   
                    break;
                end
            end
            if other_PC_break_flag==1
                break;
            end
            if last_PC_break_flag==1
                break;
            end
            cutSetSize = cutSetSize + 1;
        end
        if last_PC_break_flag==1
            break;
        end
    end
%     pc=pc_tmp;
end



    
    
    
    
    
    
    
    
    
    
% end





MB=myunion(pc,cell2mat(spouse));
% MB

time=toc(start);


