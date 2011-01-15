% NMFLAB for Signal Processing written by A. Cichocki and R. Zdunek 
% in cooperation with other members of Laboratory for Advanced Brain Signal
% Processing, BSI, RIKEN, Saitama, JAPAN

function [x,grad_x,k] = nmf_proj_grad(b,A,x,tol,no_iter,beta_dec,loss_fun,beta)

% Projected Gradient Algorithm based on the Armijo rule
%
% [x, grad_x, k]=nmf_proj_grad(b,A,x,tol,no_iter,beta_dec,loss_fun,beta) finds
% such x that solves the equation Ax = b.
%
% INPUTS:
% b - data matrix of dimension [m by T]
% A - fixed matrix of dimension [m by r]
% x - initial guess
% tol - tolerance
% no_iter - maximum number of iterations
% beta_dec - multiplier of learning rate in descent gradient direction
% loss_fun - index for the loss function (1--4)
% beta - parameter "alpha" in the Amari alpha-divergence
% 
% OUTPUTS:
% x - estimated matrix of dimension [r by T]
% grad_x - gradient of x
% k - number of performed iterations 
%
% #########################################################################

if beta == 0
   loss_fun = 2;
end
   
sigma = 0.01; alpha = 1; kinc = 1;

for k=1:no_iter  
   
 switch loss_fun  
      
      case 1 % Frobenius
          
           res = b - A*x;   
           grad_x = -A'*res;
           F = .5*norm(res,'fro')^2;
           lower_bound = 0;
                  
      case 2 % KL
                        
           bx = A*x; 
           grad_x = A'*(1 - b./(bx + eps));
           F = sum(sum(b.*log(b./(bx + eps)) + bx - b));
           lower_bound = 0.1;
           
      case 3 % Dual KL
          
           bx = A*x; 
           grad_x = A'*log(bx./b);
           F = sum(sum(bx.*log(bx./b) + b - bx));
           lower_bound = 1E-4;
                      
      case 4 % Amari alpha divergence
          
           bx = A*x; 
           grad_x = (1/beta)*A'*(1 - (b./bx).^beta);
           F = sum(sum(b.*((b./bx).^(beta - 1) - 1)/(beta^2 - beta) + (bx - b)/beta));
           lower_bound = 1E-2;
           
  end  % switch 
   
 
  tol_grad = norm(grad_x(grad_x < lower_bound | x > lower_bound));
  if tol_grad < tol,
     break
  end

%  alpha = 1;
  s = 0;
  
% search step size 
while (alpha > 1E-15) & (alpha < 1E8)
  
  s = s + 1; 
 
  xn = max(x - alpha*grad_x,lower_bound); 
  delta = xn-x; 
  
   switch loss_fun  

      case 1 % Frobenius
          
           Fn = .5*norm(b - A*xn,'fro')^2;
                  
      case 2 % KL
           
           bx = A*xn; 
           Fn = sum(sum(b.*log(b./(bx + eps)) + bx - b));
            
      case 3 % Dual KL
            
           bx = A*xn; 
           Fn = sum(sum(bx.*log(bx./b) + b - bx));
          
      case 4 % Amari alpha divergence
          
           bx = A*xn; 
           Fn = sum(sum(b.*((b./bx).^(beta - 1) - 1)/(beta^2 - beta) + (bx - b)/beta));
          
  end  % switch 
  
  cond_alpha = ((Fn - F) <= sigma*grad_x(:)'*delta(:));
  
  if cond_alpha | (xn == x)
   % Increase
     if ~kinc & (s > 1), x = xn; break; end
     alpha = alpha/beta_dec;
     kinc = 1; xp = xn;
  else
   % Decrease
     if kinc & (s > 1), x = xp; alpha = alpha*beta_dec; break; end
     alpha = alpha*beta_dec;
     kinc = 0;
  end
  
  
end % while
  
end % for k
% 
