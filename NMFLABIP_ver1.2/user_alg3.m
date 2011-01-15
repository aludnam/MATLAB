%THIS FUMCTION ILLUSTRATES HOW TO USE THE USER-DEFINED ALGORITHM
function [AH,XH] = user_alg3(Y,r,X)
%
% The example of implementing the user-defined algorithm (this is the Lee-Seung algorithm based on the Frobenius norm)
%
% INPUTS:
% Y - mixed signals (matrix of size [m by T])
% r - number of estimated signals
% X - true source signals
%
% OUTPUTS
% AH - estimated mixing matrix (matrix of size [m by r])
% XH - estimated source signals (matrix of size [r by T])
%
% #########################################################################
% Initialization
[m,T]=size(Y);
Y(Y <=0) = eps; % this enforces the positive value in the data  
AH=rand(m,r);
XH=rand(r,T);

IterNo = 1000; % number of alternating steps

% Iterations
for k = 1:IterNo             
    XH = XH.*((AH'*Y)./((AH'*AH)*XH + eps));
    AH = AH.*((XH*Y')./((AH*(XH*XH'))' + eps))';
    AH = AH*diag(1./(sum(AH,1) + eps));
end
  
            
