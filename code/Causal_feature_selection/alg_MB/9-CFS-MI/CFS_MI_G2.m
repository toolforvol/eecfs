
function   [mb,ntest,time]=CFS_MI_G2(Data,target,alpha,ns,p,k)

start=tic;
ntest=0;
mb=[];
% time=0;
% pc = yingshe2(pc,target)
%ntest=ntest+ntest1;
train_data=Data(:,mysetdiff(1:p,target));
featureMatrix = train_data;
train_label=Data(:,target);
classColumn = train_label;
numFeatures = size(featureMatrix,2);
% classScore
classScore = zeros(numFeatures,1);
vis = zeros(p,1);
for i = 1:numFeatures
    ntest=ntest+1;
%     iXY(i) = mi(featureMatrix(:,i),classColumn);
    classScore(i) = SU(featureMatrix(:,i),classColumn);
end
[classScore, indexScore] = sort(classScore,1,'descend');
% [iXY, iXYzhi] = sort(iXY,1,'descend');

% th2 = 0.0001;
threshold = 0.05;


th3 = 0.15;
t = indexScore(classScore < threshold);
u = classScore(classScore < threshold);
indexScore = indexScore(classScore > threshold);
classScore = classScore(classScore > threshold);
if ~isempty(indexScore)
    curPosition = 1;
else
    curPosition = 0;
    selectedFeatures=[];
    time=toc(start);
    return;
end

mii = -1;
while curPosition <= length(indexScore)
    mb_tmp = [];
    j = curPosition + 1;
    curFeature = indexScore(curPosition);
    while j <= length(indexScore)
        
        scoreij = SU(featureMatrix(:,curFeature),featureMatrix(:,indexScore(j)));
        ntest=ntest+1;
        if scoreij > classScore(j)
            indexScore(j) = [];
            classScore(j) = [];
%             mb_tmp = myunion(mb_tmp,indexScore(j));
        else
%             mii = classScore(j);
            j = j + 1;
        end
    end
    curPosition = curPosition + 1;
end
selectedFeatures = indexScore;
pc = selectedFeatures;
last = selectedFeatures(end);
mb = yingshe2(pc,target);
for i =1:length(selectedFeatures)
    vis(selectedFeatures(i))=1;
end
% ≈‰≈º



len1 = length(selectedFeatures);
for i=1:len1
    mb_tmp=[];
    a = find(selectedFeatures == last)+1;
    len2 = length(t);
    while a <= len2
        
        if vis(t(a))==1%ismember(t(a),mb_tmp)
            a = a + 1;
            continue;
        end
%         if ~vis(selectedFeatures(i),t(a))
%             vis(selectedFeatures(i),t(a)) = SU(featureMatrix(:,selectedFeatures(i)),featureMatrix(:,t(a)));
%             vis(t(a),selectedFeatures(i)) = vis(selectedFeatures(i),t(a));
%         end
        scoreij = SU(featureMatrix(:,selectedFeatures(i)),featureMatrix(:,t(a)));   
        ntest=ntest+1;
        if scoreij > u(a)+ 0.13
%             mb_tmp = myunion(mb_tmp,t(a));
            iXYZ = cmi(train_label,featureMatrix(:,t(a)),featureMatrix(:,selectedFeatures(i)));
            iXY0 = mi(train_label,featureMatrix(:,t(a)));
            ntest=ntest+2;
            if t(a)>=target
                ttt = t(a)+1;
            else
                ttt = t(a);
            end 
            if iXYZ > iXY0
                mb_tmp = myunion(mb_tmp,t(a));
                mb = myunion(mb,ttt);
                vis(t(a))=1;
            end    
        end
        a = a + 1;
    end
end


time=toc(start);

end

function [score] = SU(firstVector,secondVector)
%function [score] = SU(firstVector,secondVector)
%calculates SU = 2 * (I(X;Y)/(H(X) + H(Y)))
hX = h(firstVector);
hY = h(secondVector);
iXY = mi(firstVector,secondVector);

score = (2 * iXY) / (hX + hY);
end


function pc = yingshe2( pc ,target)
for i=1:length(pc)
    if pc(i)>=target
        pc(i)= pc(i)+1;
    end
end

end

        






