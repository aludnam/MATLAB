%THIS FUMCTION ILLUSTRATES HOW TO USE THE USER-DEFINED ALGORITHM
%
% This is the second-order NMF algorithm based on the Newton method. 
% The algorithm is presented by R. Zdunek and A. Cichocki in Signal
% Processing (2007)
function [A,X] = user_alg1(Y,R,Xexact)
%
% The example of implementing the user-defined algorithm 
% (this is the second-order algorithm applied to the Euclidean function)
%
% INPUTS:
% Y - mixed signals (matrix of size [M by T])
% R - number of estimated signals
% Xexact - true source signals
%
% OUTPUTS
% A - estimated mixing matrix (matrix of size [M by R])
% X - estimated source signals (matrix of size [R by T])
%
% #########################################################################
% Initialization
MaxIter = 1000; % Number of alternating steps
Alpha0 = 100; Tau = 50; % Parameters in regularization for estimation of X
%
Y(Y <=0) = eps; % this enforces the positive value in the data  
% Scaling to unit-variance
mx = 1./sqrt(var(Y,0,2) + eps);  Y = repmat(mx,[1,size(Y,2),1]).*Y;

% Settings
[M,T]=size(Y); I = eye(R); A = []; X = [];
lambda = 1E-12; % Levenberg-Marquardt regularization of Hessian
HA = spalloc(M*R,M*R,M*R^2); % Spase allocation for Hessian
while 1 A = rand(M,R); if cond(A) < 50 break; end; end % Initialization

% Alternatings
k = 0;
while k <=MaxIter
    k = k + 1;
    
  % Estimation of X  
    alpha_reg = Alpha0*exp(-k/Tau); % exponential rule for regular.param.
    if isnan(A)  disp('Matrix A is too much ill-conditioned. Try again.');
        break; end
    if cond(A) > 1E6 alphaX = 1E-6; else alphaX = 0; end
    X = max(1E6*eps,pinv(A'*A + alpha_reg + alphaX*I)*A'*Y); %Updating of X    
    
  % Estimation of A          
    hA = X*X';
    hA = hA + (lambda + 1E8*exp(-k))*eye(R);  % Levenberg-Marquardt regularization
    GA = X*Y' - hA*A';  % Gradient
    HA = kron(speye(M),-hA); % Hessian
    
cond_HA = condest(HA);
if cond_HA > 1E17
fprintf(1,'Hessian is singular in %d iteration(s). Restart is needed.\n',k) 
break; end
A = A'; A(:) = A(:) - .9*inv(HA)*GA(:); % Newton iterations
A(A <= 0) = 1E2*eps;  A = A';
A = A*diag(1./sum(A,1)); % Normalization to unit L1-column norm
end
A = repmat(1./mx,[1,size(A,2),1]).*A; % De-scaling   
    
  
            
