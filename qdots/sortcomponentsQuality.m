function [sx, isx, fsx] = sortcomponentsQuality(x,peval,rnet)
% [sx, isx, fx] = sortcomponentsQuality(x, peval, net)
% Sorts the inut values x (column vectores - like res.w) according to the
% their quality. Uses net form glm linear classifier (classified into calasses 0(rubbish),1(good),2(two good),3(three good),4(multiple good),5(half good). 
% sx: sorted compontnts;
% isx: index of sorting - x(:,isx) are sorted...
% fsx: sorted objective function used for sorting

features = computeFeatures(x, peval);
featuresW = whitenCol(features);
fPredTest=glmfwd(rnet,featuresW); % posterior densities for predicted classes
[values,labelPredTestPlusOne]=max(fPredTest,[],2);
labelPredTestOld=labelPredTestPlusOne-1;
% order labels as: 1 5 2 3 4 0 -> 1 2 3 4 5 6
labelPredTest=ones(size(labelPredTestOld));
labelPredTest(labelPredTestOld==5)=2;
labelPredTest(labelPredTestOld==2)=3;
labelPredTest(labelPredTestOld==3)=4;
labelPredTest(labelPredTestOld==4)=5;
labelPredTest(labelPredTestOld==0)=6;
% join labels with their posterior values, because then it can be sorted
% such that label 1 with highest posterior will be first, label 1 with
% second highest posterior will be second and label zero with lowest
% posterior will be last <- maybe this is not too smart...
fx = labelPredTest+(1-values); 
[sx, isx] = sort(fx, 'ascend');
fsx=sort(fx,'ascend');
