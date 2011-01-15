% GPCG algorithm adapted for NMF
function X = gpcg_nmf(Y,A,X,MaxIter,CostFun,Alpha)
% Our implementation of GPCG;

[R,T] = size(X); [M,T] = size(Y);
beta = 0.5; mu = 0.01; MaxSearchIter = 10;
Diff_DF_max = 0; gamma_GP = 0.01;

if CostFun == 1
   H = spalloc(R*T,R*T,T*R^2);
   Hx = A'*A;
   H = kron(speye(T),Hx); % Hessian
end

for k = 1:MaxIter          
    
% Step 1
    Xp = X;
    if CostFun ~= 1
       H = CostFunHess(Y,A,X,Alpha); % Hessian
    end
    P = -CostFunGrad(Y,A,X,CostFun,Alpha); % Gradient Matrix
    p = P(:); % Vectorization
    
% Step 2    
     eta = (p'*p)/(p'*(H*p));   
     for i = 0:MaxSearchIter % Armijo rule
         eta = eta*beta.^i;
         X = max(X + eta*P,0);
         if (CostFunEval(Y,A,X,CostFun,Alpha) - CostFunEval(Y,A,Xp,CostFun,Alpha)) <=  (norm(X - Xp,'fro')*mu/eta)
             break;
         end
     end
  
% Step 3     
  Z = zeros(size(X));
  Z(X > eps) = 1;
  z = vec(Z);

% Step 4  
  Pc = CostFunGrad(Y,A,X,CostFun,Alpha); % Gradient Matrix
  pc = Pc(:); % Vectorization
  
% Step 5  
  pR = z.*pc; % Reduced gradient
  if CostFun ~= 1
     H = CostFunHess(Y,A,X,Alpha); % Hessian
  end
   
HR = repmat(z,1,R*T).*H.*repmat(z',R*T,1) + speye(R*T) - spdiags(z,0,R*T,R*T); % reduced Hessian
 
% Step 6
 [pc,flag,relres,iter,resvec] = pcg(HR,-pR);
 P = reshape(pc,R,T); % Matricization
 
 % Step 7    
     DF = CostFunEval(Y,A,X,CostFun,Alpha);
     eta = 1;   
     for j = 0:MaxSearchIter % Armijo rule
         eta = eta*beta.^j;
         X = max(X + eta*P,0);
         if (CostFunEval(Y,A,X,CostFun,Alpha) < DF) 
             break;
         end
     end
     
   Fn(k) = CostFunEval(Y,A,X,CostFun,Alpha);  
     
 % Stopping criterion
   DF_old = norm(Y - A*Xp,'fro');  DF = norm(Y - A*X,'fro');
   Diff_DF = DF_old - DF;  Diff_DF_max = max(Diff_DF,Diff_DF_max);
   if (Diff_DF <= gamma_GP*Diff_DF_max) & (k > 1)
       break;
   end
        
end

% Cost Function 
function F = CostFunEval(Y,A,X,CostFun,Alpha);

 switch CostFun
     
     case 1 % Eucliden distance
         
         F = norm(Y - A*X,'fro');  
              
     case 2 % Alpha divergence
         
         Z = A*X + 1E2*eps;
         Y = Y + 1E2*eps;
         if Alpha == 1 % KL divergence
            F = sum(sum(Y.*log(Y./Z) + Z - Y));   
         elseif Alpha == 0 % Dual KL divergence
            F = sum(sum(Z.*log(Z./(Y+eps)) + Y - Z));    
         else % Alpha-divergence
            F = (1/Alpha)*sum(sum(Y.*( (Y./Z).^(Alpha - 1) - 1)/(Alpha - 1) + Z - Y));
         end
 end
   
% Gradient of Cost Function 
function G = CostFunGrad(Y,A,X,CostFun,Alpha);

 switch CostFun
     
     case 1 % Eucliden distance
         
         G = A'*(A*X - Y);
              
     case 2 % Alpha divergence
         
         Z = A*X+1E2*eps;
         if ~Alpha % Dual KL divergence
            G = A'*log(Z./(Y + 1E2*eps));
         else
            G = (1/Alpha)*A'*(1 - ((Y+1E2*eps)./Z).^Alpha);
         end
  end

% Hessian of Cost Function 
function H = CostFunHess(Y,A,X,Alpha);

[R,T] = size(X);
M = size(Y,1);
H = spalloc(R*T,R*T,T*R^2);
Z = A*X+1E2*eps;
if ~Alpha % Dual KL divergence
    Zx = 1./Z + 1E2*eps;
else
    Zx = ((Y+1E2*eps).^Alpha)./(Z.^(Alpha + 1));
end
 for t = 1:T
     H(((t-1)*R+1):t*R,((t-1)*R+1):t*R) = A'*repmat(Zx(:,t),1,M)*A;
 end
        
            
            
            
            
            