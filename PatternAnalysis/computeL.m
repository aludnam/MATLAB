function L = computeL(xK,K)
% COMPUTEL computes L-function from ripleys K function (L=sqrt(K/pi)-xK)
% L = computeL(xK,K)
% xK - position where K function was estimated
% K - K-function
L=bsxfun(@minus,sqrt(K/pi),xK);