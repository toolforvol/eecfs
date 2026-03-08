function [MB,test,time] = EECFS_Z(Data,target,alpha,ns,p,maxK)

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
dep_sp=zeros(p,p);
spouse = cell(1,p);

sp=[];
MB=[];

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
PCselect = myunion(PCinter, PCdiff);

% step4:find the spouse nodes(different from EDMB)

NoPC = mysetdiff(1:p,[PCselect,target]);

for i =1:length(PCselect)
    A = PCselect(i); % A equals to X in the CFS framework of the literature
    for j=1:length(NoPC)
        Y = NoPC(j); 

        if ~isempty(find(dSep{Y}==A, 1))
            continue
        end 

        test = test+1;
        [pval,dep_sp(A,Y)]=my_fisherz_test(Y,target,myunion(dSep{Y},A),Data,ns,alpha);   
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
            test = test+1;
            [pval]=my_fisherz_test(Y,A,[],Data,ns,alpha); 

            if isnan(pval)
                pval=0;
            end

            if pval<=alpha
                spouse{1,A} = myunion(spouse{1,A},Y); 
            end

        end   

    end

end


for i =1:length(PCselect)
    A = PCselect(i);

    [spouse{1,A},test1] = MMPC_Z_1224(Data,A,alpha,ns,p,maxK,spouse{1,A}); 

    test = test+test1;
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

    spouse{1,A}=tmp_PC;

end

%fid = fopen(['EECFS_RNAnew_spouse_links_target' num2str(target) '.txt'], 'w');
%fprintf(fid, 'Y\tA\tTarget\n');
%for A = 1:p
%    if ~isempty(spouse{1,A})
%        for y = spouse{1,A}
%            fprintf(fid, '%d\t%d\t%d\n', y, A, target);
%        end
%    end
%end
%fclose(fid);


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
disp(PCselect)

time=toc(start);

end


% function [MB,test,time] = MEDMB_Z(Data,target,alpha,ns,p,maxK)
% 
% %
% % MEDMB_Z finds the Markov blanket of target node on discrete data
% % (modified to log spouse discoveries similarly to MMMB_Z)
% %
% % INPUT :
% %       Data is the data matrix
% %       target is the index of target node
% %       alpha is the significance level
% %       ns is the size array of each node
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
% if (nargin == 3)
%    ns=max(Data);
%    [~,p]=size(Data);
%    maxK=3;
% end
% 
% start=tic;
% test = 0;
% dep_sp=zeros(p,p);
% spouse = cell(1,p);
% 
% sp=[];
% MB=[];
% 
% % open log file (include target in filename)
% logFileName = sprintf('MEDMB_PROTEIN_causal_structure_log_target_%d.txt', target);
% logFile = fopen(logFileName, 'w');
% if logFile ~= -1
%     fprintf(logFile, '=== Causal Feature Selection Log (MEDMB_Z) ===\n');
%     fprintf(logFile, 'Target Node: %d\n\n', target);
% else
%     warning('Could not open log file %s for writing. Proceeding without logging to file.', logFileName);
% end
% 
% % step1:find the candidate PC sets
% % Use PCMB and MBOR algorithm to obtain PC sets
% 
% [PCMB_CPC,test1,~,dSep] = PCMB_PC_Z(Data,target,alpha,ns,p,maxK);
% [MBOR_CPC,test2,~] = MBOR_PC_Z(Data,target,alpha,ns,p,maxK);
% 
% test = test1 + test2;
% 
% if logFile ~= -1
%     fprintf(logFile, '[Step1] PCMB_CPC of target %d: %s\n', target, mat2str(PCMB_CPC));
%     fprintf(logFile, '[Step1] MBOR_CPC of target %d: %s\n\n', target, mat2str(MBOR_CPC));
% end
% 
% % step2:obtain the intersection,union and difference of the PC sets
% PCinter = myintersect(PCMB_CPC,MBOR_CPC);
% PCunion = myunion(PCMB_CPC,MBOR_CPC);
% PCdiff = mysetdiff(PCunion,PCinter);
% 
% if logFile ~= -1
%     fprintf(logFile, '[Step2] PCinter: %s\n', mat2str(PCinter));
%     fprintf(logFile, '[Step2] PCunion: %s\n', mat2str(PCunion));
%     fprintf(logFile, '[Step2] PCdiff (before removal): %s\n\n', mat2str(PCdiff));
% end
% 
% % step3:remove wrong PC nodes from PCdiff
% 
% %NonTPC euqals to U\PCunion\{T}
% NonTPC = mysetdiff(1:p,myunion(PCunion,target)); 
% remove = [];
% 
% for i = 1:length(PCunion)
%     Y = PCunion(i);
%     break_flag=0;
%     for j=1:length(NonTPC)
%         X = NonTPC(j);
%         % S is the conditoned set or separated set
%         S =  myunion(dSep{1,X}, Y);   
%         test = test +1;
%         % Judge independence or conditional independence among variables
%         [pval]=my_fisherz_test(target,X,S,Data,ns,alpha); 
%         if isnan(pval)        
%             CI=0;
%         else
%             if pval<=alpha
%                 CI=0;
%             else
%                 CI=1;
%             end
%         end
%         % According to independence result,remove wrong PC nodes from PCdiff
%         if(CI==0)    
%            SZ = mysetdiff( myunion(PCunion,X),Y );
%            cutSetSize = 0;
%            while length(SZ) >= cutSetSize&&cutSetSize<=maxK
%                SS = subsets1(SZ, cutSetSize);   
%                for si=1:length(SS)
%                    S = SS{si};
%                    test=test+1;
%                    [pval]=my_fisherz_test(Y,target,S,Data,ns,alpha); 
%                    if isnan(pval)      
%                        CI=0;
%                    else
%                        if pval<=alpha
%                            CI=0;
%                        else
%                            CI=1;
%                        end
%                    end
%                    if(CI==1)  
%                        remove = myunion(remove,Y);
%                        break_flag=1; 
%                        break;
%                    end
%                end
%                if( break_flag==1 )
%                    break;
%                end
%                cutSetSize = cutSetSize + 1;       
%            end
%            if( break_flag==1 )
%                break;
%            end
%         end
%     end
% end
% % Process PCdiff set
% PCdiff = mysetdiff(PCdiff,remove);
% 
% % Obtain PCselect
% PCselect = myunion(PCinter, PCdiff);
% 
% if logFile ~= -1
%     fprintf(logFile, '[Step3] PCselect: %s\n\n', mat2str(PCselect));
% end
% 
% % step4:find the spouse nodes(different from EDMB)
% 
% NoPC = mysetdiff(1:p,[PCselect,target]);
% 
% for i =1:length(PCselect)
%     % A also equals to X in the CFS framework
%     A = PCselect(i); % A equals to X in the CFS framework of the literature
%     for j=1:length(NoPC)
%         Y = NoPC(j); %Y is current node not in PC of target
% 
%         if ~isempty(find(dSep{Y}==A, 1))
%             continue
%         end 
% 
%         test = test+1;
%         [pval,dep_sp(A,Y)]=my_fisherz_test(Y,target,myunion(dSep{Y},A),Data,ns,alpha);   
%         if isnan(pval)
%             CI=1;
%         else
%             if pval<=alpha
%                 CI=0;
%             else
%                 CI=1;
%             end
%         end
% 
%         if CI==0 
%             test = test+1;
%             [pval]=my_fisherz_test(Y,A,[],Data,ns,alpha); 
% 
%             if isnan(pval)
%                 pval=0;
%             end
% 
%             if pval<=alpha
%                 spouse{1,A} = myunion(spouse{1,A},Y); %
%                 if logFile ~= -1
%                     fprintf(logFile, '[Spouse-Initial] Potential spouse found: Y=%d is potential spouse of Target=%d via A=%d (added to spouse{%d}).\n', Y, target, A, A);
%                 end
%             end
% 
%         end   
% 
%     end
% 
% end
% 
% % Second-stage refine spouse sets by learning PC of A (MMPC_Z_1224)
% for i =1:length(PCselect)
%     A = PCselect(i);
% 
%     [spouse{1,A},test1] = MMPC_Z_1224(Data,A,alpha,ns,p,maxK,spouse{1,A}); 
% 
%     test = test+test1;
%     %--------------------------------------2---------------
%     tmp_PC = spouse{1,A};
%     for j=(length(spouse{1,A})):-1:1
%         break_flag=0;
%         X = spouse{1,A}(j);
% 
%         CanPC=mysetdiff(myunion(tmp_PC,PCselect), [X,A]);
% 
%         cutSetSize = 1;
% 
%         while length(CanPC) >= cutSetSize&&cutSetSize+1<=maxK
% 
%             SS = subsets1(CanPC, cutSetSize);    
% 
%             for si=1:length(SS)
%                 Z = SS{si};
% 
%                 test=test+1;
% 
%                 [pval]=my_fisherz_test(target,X,myunion(Z,A),Data,ns,alpha);    
% 
%                 if isnan(pval)
%                     CI=0;
%                 else
%                     if pval<=alpha
%                         CI=0;
%                     else
%                         CI=1;
%                     end
%                 end
% 
%                 if (CI==1)     
%                       tmp_PC = mysetdiff(tmp_PC,X);
%                       if logFile ~= -1
%                           fprintf(logFile, '[Spouse-Refine] Removed node %d from spouse{%d} after conditioning on Z=%s (A=%d, target=%d).\n', X, A, mat2str(Z), A, target);
%                       end
%                       break_flag=1;
%                       break;
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
%     spouse{1,A}=tmp_PC;
%     if logFile ~= -1
%         fprintf(logFile, '[Spouse-Final] Final spouse{%d}: %s\n', A, mat2str(spouse{1,A}));
%     end
% 
% end
% 
% 
% % write all final spouse links to log file
% if logFile ~= -1
%     fprintf(logFile, '\n=== Final Results for target %d ===\n', target);
%     fprintf(logFile, 'PCselect: %s\n', mat2str(PCselect));
%     fprintf(logFile, 'Spouse links (Y\tA\tTarget):\n');
%     for A = 1:p
%         if ~isempty(spouse{1,A})
%             for y = spouse{1,A}
%                 fprintf(logFile, '%d\t%d\t%d\n', y, A, target);
%             end
%         end
%     end
%     fprintf(logFile, '\nMB (union of PCselect and spouses):\n');
% end
% 
% MB=myunion(PCselect,cell2mat(spouse));
% disp(PCselect)
% 
% if logFile ~= -1
%     fprintf(logFile, 'MB: %s\n', mat2str(MB));
%     fclose(logFile);
% end
% 
% time=toc(start);
% 
% end

