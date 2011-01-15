%THIS FUMCTION ILLUSTRATES HOW TO USE THE USER-DEFINED ALGORITHM
function [AH,XH] = user_alg1(Y,r,X)
%
% The example of implementing the user-defined algorithm (this is the regularized Fixed-Point algorithm based)
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
    
    alpha_reg = 20*exp(-k/10); % regularization parameter
    XH = max(1E6*eps,pinv(AH'*AH +  alpha_reg)*AH'*Y);   
    AH = max(1E6*eps, Y*XH'*pinv(XH*XH' + alpha_reg));  
    AH = AH*diag(1./(sum(AH,1) + eps));
    
end
  
            
