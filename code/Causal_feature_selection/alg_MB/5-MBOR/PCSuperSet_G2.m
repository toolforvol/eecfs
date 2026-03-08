function [PCS,test,sepset,time] = PCSuperSet_G2(Data,target,alpha,ns,p,maxK)

start=tic;
test = 0;

sepset=cell(1,p);

% ------------------------------------------------------------
% Phase I: Remove X if T ⊥ X

PCS = mysetdiff(1:p,target);
tmp_PCS = PCS;
for i=1:length(PCS)
    X = PCS(i);
    test=test+1;
    [pval]=my_g2_test(X,target,[],Data,ns,alpha);
    if pval>alpha
        tmp_PCS = mysetdiff(tmp_PCS,X);
        sepset{X} = {};%若要使用MBOR的dSep将{}改成[]即可
    end
end
PCS = tmp_PCS;

% ------------------------------------------------------------
% Phase II:Remove X if T ⊥ X|Y

tmp_PCS = PCS;
for i=1:length(PCS)
    X = PCS(i);
    PCSX = mysetdiff( PCS,X );
    for j=1:length(PCSX)
        Y = PCSX(j);
        test=test+1;
        [pval]=my_g2_test(X,target,Y,Data,ns,alpha);
        if pval>alpha
            tmp_PCS = mysetdiff(tmp_PCS,X);
            sepset{X} = Y;
            
            break;     
        end
    end
end
PCS = tmp_PCS;
time=toc(start);





