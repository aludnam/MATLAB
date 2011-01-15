function X = ipg(A,Y,X,maxiter)
%
% Interior-Point Method
%
% REFERENCE:
%   M. Merritt and Y. Zhang. An Interior-Point Gradient Method for
%   Large-Scale Totally Nonnegative Least Squares Problems. J. Optimization
%   Theory and Applications, Vol. 126. No. 1, pp. 191-202, July, 2005.
%

% initialize
B = A'*A;
BX = B*X;

for iter = 1:maxiter
   
   % step computation
    R = A*X - Y;
    S = A'*R;
    P = (X./max(BX,eps)).*S; 
    BP = B*P; 

    % steplength
     alpha1 = sum(P.*S,1)./max(sum(P.*BP,1),eps);
     alpha2 = -1./min(-eps,min(-P./max(X,eps),[],1) );
     alpha = min([.99*alpha2; alpha1],[],1);
    
    % update
    X = X - P.*repmat(alpha,size(B,2),1);
    BX = BX - BP.*repmat(alpha,size(B,2),1); 
    
end


