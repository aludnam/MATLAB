function [l,fPredTest]=labelW(w,peval,rnet)
% [lablelPred,fPred]=labelW(w,peval,rnet)
% Labels data from the linear classifier specified by rnet. (or loads the
% one trained on real data.)
% w = nx*ny X ncomp
% rnet = linear classifier
% labelPred - predicted labels: 1-good,2=two good, 3-4 multiple good,5-good half missing, 0 - rubbish.
% fPred - posterior over classes

if ~exist('rnet','var')
    load ~/project/data/qdots/S371/rnet
end
features = computeFeatures(w, peval);
featuresW = whitenCol(features);
fPredTest=glmfwd(rnet,featuresW); % posterior densities for predicted classes
[values,labelPredTestPlusOne]=max(fPredTest,[],2);
l=labelPredTestPlusOne-1;