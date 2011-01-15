% Nonnegative Matrix Factorization with Quadratic Programming
% The function iptrnr(Y,X,no_iter,lambdaA,rho,eta,epsil_A) returns
% the matrix A for the constrained system of linear equations: AX = Y

function A = iptrnr(Y,X,no_iter,lambdaA,rho,eta,epsil_A)

% INPUTS:
% > Y:        observation matrix [I by K] of mixed signals (images) 
% > X:        source component matrix [J by K]
% > no_iter:  number of inner iterations
% > lambdaA:  Tikhonov regularization parameter
% > rho:      auxiliary parameter for defining the log barrier parameter
% > eta:      parameter in stopping rule for iterations
% > epsil_A:  threshold for active-set variables

% OUTPUT:
% > A:        mixing matrix [I by J]

[J,K] = size(X); [I,Ky] = size(Y); M = I*J;

% Initialization
A = ones(I,J);
Q_bar = spalloc(M,M,M*J); Q_tilde = spalloc(M,M,M*J);

B = X*X'; C_bar = X*Y';
Q_bar = kron(speye(I),B); % Hessian of D(Y||AX) with respect to A
Go = (C_bar - B*A'); % Gradient of D(Y||AX) with respect to A^T
At = A';
theta = -(rho/(M))*abs(Go(:)'*At(:)); % initialization for 
                                      %log barrier parameter
l = 0;
while l < no_iter % inner iterations
    
    l = l + 1; 
    a = At(:);
    active = find(a < epsil_A); % active set
    inactive = find(a >= epsil_A); % inactive set
     
    if ~isempty(active)
        a_tilde = a;             a_tilde(active) = [];
        N = length(a_tilde);
        c_tilde = C_bar(:);      c_tilde(active) = [];
        c_tilde = 2*theta*1./a_tilde - c_tilde;
        Q_tilde = Q_bar;
        Q_tilde(:,active) = [];  Q_tilde(active,:) = [];
        Q_tilde = Q_tilde - theta*spdiags(1./a_tilde.^2,0,N,N);
    else
        a_tilde = a;
        Q_tilde = Q_bar - theta*spdiags(1./a.^2,0,M,M) + lambdaA*speye(M);
        c_tilde = 2*theta*(1./a) - C_bar(:);
    end
   
     [Q,R] = qr(Q_tilde,-c_tilde); % Q-less QR factorization
     z_tilde = R\Q; % Gaussian elimination
     h_tilde = z_tilde - a_tilde;
     beta = min(1, 0.9995*min(a_tilde./abs(h_tilde)));
     a_tilde = a_tilde + beta*h_tilde;
     
     theta_tilde = theta*(2*(1./a_tilde) - (1./a_tilde.^2).*z_tilde);
     theta = ( theta_tilde'*a_tilde);
         
     if ~isempty(active)
        an = zeros(M,1);   an(inactive) = a_tilde;   a = an;
     else
        a = a_tilde; 
     end
    
     if abs(theta) < eta*norm(a) % stopping rule
        break
     end
     theta = (rho/sqrt(M))*theta;
     At = reshape(a,J,I);
         
end
A = At';



