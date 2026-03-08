
function [MB,test,time] = CFS_GetPCD(Data,target,alpha,ns,p,maxK)

if (nargin == 3)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end


start=tic;

spouse = cell(1,p);

dep_sp=zeros(p,p);

% Remove the nodes independent of the target node conditioning on
% empty set, and sort the dependent nodes

score=1000000;

% % 

[pc,test,~,sepset] = GetPCD_G2_1224(Data,target,alpha,ns,p,maxK,mysetdiff(1:p,target));
NoPC = mysetdiff(1:p,[pc,target]);


for i =1:length(pc)
    X = pc(i);
    % CFS框架中的5-9行
    for j=1:length(NoPC)
        Y = NoPC(j);
        if ~isempty(find(sepset{Y}==X, 1))
            continue
        end 
        
        test = test+1;
        [pval,dep_sp(X,Y)]=my_g2_test(Y,target,myunion(sepset{Y},X),Data,ns,alpha);    

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
            spouse{1,X} = myunion(spouse{1,X},Y); %现在得到的spouse{1,X}就是CFS框架中第七行的CSPT{X}
        end   
        
    end
    
    A=pc(i); 
    [spouse{1,A},test1] = GetPCD_G2_1224(Data,A,alpha,ns,p,maxK,spouse{1,A}); 
    test = test + test1;
    
    SP=[];
    [~,SP_index]=sort(dep_sp(A,spouse{1,A}),'descend'); 
    for f=1:length(spouse{1,A})
        B = spouse{1,A}(SP_index(f));
        SP=[SP B]; 
        SP_length=length(SP);
        SP_tmp=SP;
        SP_break_flag=0;
        for c=SP_length:-1:1 
            X = SP(c); 
            CanSP=mysetdiff(myunion(pc,SP_tmp), [X,A]);
            cutSetSize = 1;            
            
            other_SP_break_flag=0;
            while length(CanSP) >= cutSetSize&& cutSetSize+1<=maxK
                SS = subsets1(CanSP, cutSetSize);  
                for si=1:length(SS)
                    Z = SS{si};
                    condset=myunion(Z,A);
                    
                    if (length(intersect(pc,condset))==1)
                        continue;
                    end                    
                    
                    if B~=X    
                        if isempty(find(Z==B, 1))
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
                        X
                        condset
                        SP_tmp=mysetdiff(SP_tmp, X);
                        
                        if B==X %alter
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
    spouse{1,A} = SP;
end




% % %--------------------------------------- remove false
% % 
% % % Phase I: Remove false positives from SPST
% % 

% % % Phase II: Remove false positives from PCST
% % 
% % 

pval=ones(1,p)*score;
dep_tmp1 = ones(1,p)*score;  
last_added=-1;

CanPCD = pc;
pc=[];
dep_re = ones(1,p)*score;

dep_tmp_remove = ones(1,p)*score;
pval_remove=ones(1,p)*score;

while( ~isempty(CanPCD) )
    Sep=cell(1,p);
    
    dep=zeros(1,score);
    dep_tmp2 = -score;
    

    tmp_pval=ones(1,p)*score;
    Y=-1;
    
    %----------------GET PC STEP 1-1: get Sep[X]
    
    for i=1:length(CanPCD)
                       
        X = CanPCD(i);
        cutSetSize = 0;
        
        
        while length(pc) >= cutSetSize &&cutSetSize<=maxK
            
            SS = subsets1(pc, cutSetSize);   
            
            for si=1:length(SS)
                Z = SS{si};
                
                if last_added~=-1
                    if ~ismember(last_added,Z)
                        continue;
                    end
                end
                
                test=test+1;
                
                [tmp_pval(X),dep(si)]=my_g2_test(X,target,Z,Data,ns,alpha);     
                
                if (dep(si) < dep_tmp1(X))           
                    
                    dep_tmp1(X) = dep(si);
                    dep_re(X)=dep_tmp1(X);
                    Sep{X} = Z;
                    pval(X)=tmp_pval(X);
                end
            end

            cutSetSize = cutSetSize + 1;
        end
        
    end
    
    tmp_CanPCD = CanPCD;
    for i=1:length(CanPCD)
                       
        X = CanPCD(i);
        if isnan(pval(X))
            CI=0;
        else
            if pval(X)<=alpha
                CI=0;
            else
                CI=1;
            end
        end
        
        if CI==1
            sepset{X}=Sep{X};
            NoPC = [NoPC,X];
            tmp_CanPCD = mysetdiff(tmp_CanPCD,X);
        end
    end
    CanPCD = tmp_CanPCD;

    %----------------GET PC STEP 1-1:  END
    
    %----------------GET PC STEP 1-2:  get Y
    
    for i=1:length(CanPCD)
        X = CanPCD(i);
        if dep_tmp1(X)>dep_tmp2&&dep_tmp1(X)~=score
            Y=X;
            dep_tmp2=dep_tmp1(X);
        end
    end
    %----------------GET PC STEP 1-2:  END
        
    %----------------GET PC STEP 1-3:  test CI T with Y 

    if Y~=-1
        if pval(Y)<=alpha||isnan(pval(Y))
            pc=[pc Y];
            CanPCD = mysetdiff(CanPCD,Y);
        else
            break;
        end
    else
        break;
    end

    last_added=Y;
    
    %----------------GET PC STEP 1-3:  END
    
    
    % -----------------------------------------------------------
    % remove false positives from PCD
           
    Sep_remove=cell(1,p);
    
    dep_remove=zeros(1,p);   

    tmp_pval_remove=ones(1,p)*score;
    
    for i=1:(length(pc))

        
                
        X = pc(i);

        cutSetSize = 0;
        condset = mysetdiff(pc,X);
               
        while length(condset) >= cutSetSize &&cutSetSize<=maxK
            
            SS = subsets1(condset, cutSetSize);   
            
            for si=1:length(SS)
                Z = SS{si};

                if ~ismember(last_added,Z)&&X~=last_added
                    continue;
                end
                spouse_test=[];
                for k=1:length(Z)
                    pc_var = Z(k);
                    spouse_test = myunion(spouse_test,spouse{1,pc_var});
                end
                TestSet = myunion(Z,spouse_test);                
                
                
                test=test+1;
                [tmp_pval_remove(X),dep_remove(si)]=my_g2_test(X,target,TestSet,Data,ns,alpha);     
                

                if (dep_remove(si) < dep_tmp_remove(X))          
                    
                    dep_tmp_remove(X) = dep_remove(si);
                    dep_re(X) = dep_tmp_remove(X);
                    Sep_remove{X} = Z;
                    pval_remove(X)=tmp_pval_remove(X);
                end
                
            end
            
            cutSetSize = cutSetSize + 1;
        end
        
    end
    
    tmp_PCD = pc;
    
    for i=1:(length(pc))
                        
        X = pc(i);
        
        if pval_remove(X)>alpha %&& pval_remove(X)~=score   %alter---------
            sepset{X}=Sep_remove{X};
            tmp_PCD = mysetdiff( tmp_PCD,X );
            NoPC = [NoPC,X];
        end
    end
    
    pc = tmp_PCD;

end
    

MB=myunion(pc,cell2mat(spouse));
% MB

time=toc(start);