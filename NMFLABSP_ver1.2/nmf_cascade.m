% NMFLAB for Signal Processing written by A. Cichocki and R. Zdunek 
% in cooperation with other members of Laboratory for Advanced Brain Signal
% Processing, BSI, RIKEN, Saitama, JAPAN

function [A,X,Distance_output]=nmf_cascade(Y,r,Index_norm,A_true,S,mc_on, NoAlts,restart_mc_on,AL,Y_true,type_alg,max_restart,no_iter, alphaA1, alphaA2, alphaX, alphaSa, alphaS, alpha)
%
%
% Non-negative Matrix Factorization (NMF) with the cascade implementation
% of some NMF algorithms
%
% [A,X]=nmf_cascade(Y, r, Index_norm, A_true, S, mc_on, NoAlts, restart_mc_on, AL, Y_true, type_alg, max_restart, no_iter, 
%       alphaA1, alphaA2, alphaX, alphaSa, alphaS, alpha)
%       returns mixing matrix A of dimension [m by r],
%       and source matrix X of dimension [r by T], for the linear mixing model: AX = Y, 
%       where Y is an observation matrix [m by T]. 
% Note: > m: number of sensors,
%       > r: number of sources,
%       > T: number of samples,
% 
% INPUTS:
%       > Index_norm:    vector of 13 binary entries indicating which number
%                        divergence measures are turned on in the View Options,
%
%       > A_true:        true mixing matrix (only for synthetic data),
%       > S:             true source matrix (only for synthetic data), 
%       > mc_on:         1 - Monte Carlo analysis enabled, 0 - otherwise, 
%       > NoAlts:        number of alternating steps (only for Monte Carlo
%                        analysis and the option "Fixed Alternatings" is selected)
%       > restart_mc_on: 1 - restarts in Monte Carlo analysis are enabled,
%                        0 - otherwise, 
%       > AL:            mixing matrix estimaed from the preceeding layer, 
%       > Y_true:        the first layer mixed signals (mixtures),
%       > type_alg:      indicates the selected algorithm,  
%       > max_restart:   number of restarts, 
%       > no_iter:       number of inner iterations, 
%       > alphaA1:       regularization parameter for computation of
%                        the first mixing matrix,
%       > alphaA2:       regularization parameter for computation of
%                        the second mixing matrix,
%       > alphaX:        regularization parameter for computation of
%                        the sources,
%       > alphaSa:       parameter of non-linear projection in computation
%                        of the mixing matrix,
%       > alphaS:        parameter of non-linear projection in computation
%                        of the sources,
%       > alpha:         parameter "alpha" in the Amari alpha-divergence,
% 
% OUTPUTS:
%       > A:               estimated mixing matrix,
%       > X:               estimated source matrix,
%       > Distance_output: structures of different divergences measured between "Y" and estimated "AX" versus iterations,
%
%
% #########################################################################
A = [];
X = [];
if (nargin < 19) | isempty(alpha) | max(size(alpha) > 1)
   disp('Incorrect parameter alpha in alpha-divergence');
   return
end
if (nargin < 18) | isempty(alphaS) | max(size(alphaS) > 1)
   disp('Incorrect regularization parameter alphaS');
   return
end
if (nargin < 17) | isempty(alphaSa) | max(size(alphaSa) > 1)
   disp('Incorrect regularization parameter alphaSa');
   return
end
if (nargin < 16) | isempty(alphaX) | max(size(alphaX) > 1)
   disp('Incorrect regularization parameter alphaX');
   return
end
if (nargin < 15) | isempty(alphaA2) | max(size(alphaA2) > 1)
   disp('Incorrect regularization parameter alphaA2');
   return
end
if (nargin < 14) | isempty(alphaA1) | max(size(alphaA1) > 1)
   disp('Incorrect regularization parameter alphaA1');
   return
end
if (nargin < 13) | isempty(no_iter) | (no_iter < 1) | max(size(no_iter) > 1)
   disp('Incorrect number of iterations in the EMML algorithm');
   return
end
if (nargin < 12) | isempty(max_restart)  | (max_restart < 0) | max(size(max_restart) > 1)
   disp('Number of restarts must be given correctly');
   return
end
if (nargin < 11) | isempty(type_alg) | (type_alg < 0) | max(size(type_alg) > 1)
   disp('Unknown algorithm');
   return
end
if (nargin < 10) | isempty(Y_true) 
   disp('The first layer mixed signals are unknown');
   Y_true = zeros(size(Y_true));
end
if (nargin < 9) | isempty(AL) 
   disp('Mixing matrix from the preceeding layer unknown');
   AL = eye(size(Y,1));
end
if (nargin < 8) | isempty(restart_mc_on) | max(size(restart_mc_on) > 1)
   disp('Index od restarts in MC analysis unknown');
   restart_mc_on = 0;
end
if (nargin < 7) | isempty(NoAlts) | max(size(NoAlts) > 1)
   disp('Adjustable number of alternatings');
   NoAlts = [];
end
if (nargin < 6) | isempty(mc_on) | max(size(mc_on) > 1)
   disp('No Monte Carlo Analysis');
   mc_on = 0;
end
if (nargin < 5) | isempty(S) 
   disp('X_true not given');
end
if (nargin < 4) | isempty(A_true) 
   disp('A_true not given');
   index_fixed_A = 1;
else
   index_fixed_A = 0;  
end
if (nargin < 3) | isempty(Index_norm)
   '"Index_norm" must be specified'
   return
end
if (nargin < 2) | isempty(r)
   'Rank of factorization must be given'
   return
end
if isempty(Y) | isnan(Y)
   error('No data');
   return
end
% test for negative values in Y
if min(min(Y)) < 0
    disp('Some matrix entries are changed from negative to small positive');
    Y(Y< 0) = eps;
end
if min(sum(Y,2)) == 0
    disp('Not all entries in a row can be zero');
    return
end

Y = Y + eps;
%Y = Y*Y';

if (alpha == 0) & (type_alg == 3)
    type_alg = 4;
end
   
[m,T]=size(Y);
niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)
niter_sample = 30; % maximum number of iterations for each random sample
epsil_normA = 1E-12; % tolerance for alternating

% Monte Carlo and alternatings adjustment
if mc_on & ~restart_mc_on
    max_restart = 0;
end
if ~isempty(NoAlts)
    niter_selected = NoAlts;
end

% Declaration for A and X
A=zeros(m,r);
Ap = zeros(m); 
X=zeros(r,T);
Ainit = A;
Xinit = X;
Z = zeros(m,T);
KL_outer_temp = 0;
Z_outer = 0;
nr = 0; restart_on = 0; norm_A = 10; nr_best = -1;
m_sx = 1:m; r_sx = 1:r; T_sx = 1:T; s_dist = 0;


while (nr <= max_restart)
        
   % Initialize random A and X
       if ~nr & (~mc_on | restart_mc_on) 
          A1_init(m_sx',m_sx) = abs(repmat(.1*sin(2*pi*.1*m_sx'),1,m) + repmat(.1*cos(2*pi*.1*m_sx),m,1) + repmat(cos(2*pi*.471*m_sx'),1,m) + repmat(sin(2*pi*.471*m_sx),m,1));
          A1_init = A1_init/max(max(A1_init));
        
          A2_init(m_sx',r_sx) = abs(repmat(.1*sin(2*pi*.1*m_sx'),1,r) + repmat(.1*cos(2*pi*.1*r_sx),m,1) + repmat(cos(2*pi*.471*m_sx'),1,r) + repmat(sin(2*pi*.471*r_sx),m,1));
          A2_init = A2_init/max(max(A2_init));
        
          Xinit(r_sx',T_sx) = abs(repmat(.1*sin(2*pi*.1*r_sx'),1,T) + repmat(.1*cos(2*pi*.1*T_sx),r,1) + repmat(cos(2*pi*.471*r_sx'),1,T) + repmat(sin(2*pi*.471*T_sx),r,1));
          Xinit = Xinit/max(max(Xinit));
       else
          A1_init=rand(m,m);
          A2_init=rand(m,r);
          Xinit=rand(r,T);
       end
       
        % Normalization of initial guess
       A1_init = A1_init*diag(1./sum(A1_init,1));
       A2_init = A2_init*diag(1./sum(A2_init,1));
        
        if (nr == max_restart)&(max_restart > 0)
           A1 = A1_best;
           A2 = A2_best;
           X = X_best;
        else
           A1 = A1_init;
           A2 = A2_init;
           X = Xinit;
        end % initial guess assignment
    
    Yx = zeros(m,T);
    n = 0; k = 0;
    
while ((k <= niter_sample)&(nr < max_restart)) | ((k <= niter_selected)&(nr == max_restart)&(norm_A > epsil_normA)& isempty(NoAlts)) | ((k <= niter_selected)&(nr == max_restart)& (NoAlts > 0)) 
 
k = k + 1;
    
% generalized divergence-reducing NMF iterations (main algorithm)
        if no_iter == 1
            
            if type_alg == 1 % KL
                
                Ap = A1;
               
                Nom = (Y./(A1*A2*X + eps))*(A2*X)';
                Nom(Nom <= 0) = eps;
                A1 = (A1.*(Nom./( repmat((A2*X*ones(T,1))',m,1) + eps))).^(1 + alphaSa);
                A1 = A1*diag(1./sum(A1,1));
                
                Nom = A1'*(Y./(A1*A2*X + eps))*X';
                Nom(Nom <= 0) = eps;
                A2 = (A2.*(Nom./(A1'*ones(m,1)*ones(1,T)*X' + eps))).^(1 + alphaSa);
                A2 = A2*diag(1./sum(A2,1));
                              
                Nom = (A1*A2)'*(Y./(A1*A2*X + eps));
                Nom(Nom <= 0) = eps;
                X = (X.*(Nom./(repmat((A1*A2)'*ones(m,1),1,T) + eps))).^(1 + alphaS);


            elseif type_alg == 2 % Frobenius
                
                Ap = A1; 
                Nom = Y*(A2*X)' - alphaA1;
                Nom(Nom <= 0) = eps;
                A1 = A1.*(Nom./(A1*A2*X*(A2*X)' + eps));
                A1 = A1*diag(1./sum(A1,1));
                  
                Nom = A1'*Y*X' - alphaA2;
                Nom(Nom <= 0) = eps;
                A2 = A2.*(Nom./(A1'*A1*A2*X*X' + eps));
                A2 = A2*diag(1./sum(A2,1));
                  
                Nom = (A1*A2)'*Y - alphaX;
                Nom(Nom <= 0) = eps;
                X = X.*(Nom./((A1*A2)'*A1*A2*X + eps));
                         
            elseif type_alg == 3 % Amari alpha-divergence
                
                Ap = A1; 
                Nom = ((Y./(A1*A2*X + eps)).^alpha )*(A2*X)';
                Nom(Nom <= 0) = eps;
                A1 = (A1.*(Nom./( repmat((A2*X*ones(T,1))',m,1) + eps)).^(1/alpha) ).^(1 + alphaSa);
                A1 = A1*diag(1./sum(A1,1));
                
                Nom = A1'*((Y./(A1*A2*X + eps)).^alpha )*X';
                Nom(Nom <= 0) = eps;
                A2 = (A2.*(Nom./(A1'*ones(m,1)*ones(1,T)*X' + eps)).^(1/alpha) ).^(1 + alphaSa);
                A2 = A2*diag(1./sum(A2,1));
                
                Nom = (A1*A2)'*(Y./(A1*A2*X + eps)).^alpha;
                Nom(Nom <= 0) = eps;
                X = (X.*(Nom./(repmat((A1*A2)'*ones(m,1),1,T) + eps)).^(1/alpha) ).^(1 + alphaS);
            
            elseif type_alg == 4 % SMART
                
                Ap = A1; 
                Nom = ( log(Y./(A1*A2*X + eps) +eps) )*(A2*X)';
                A1 = (A1.*exp( Nom* diag(1./(sum( A2*X,2) +eps) ) ) ).^(1 + alphaSa);
                A1 = A1*diag(1./(sum(A1,1) + eps));
                
                Nom = A1'*( log(Y./(A1*A2*X + eps) + eps) )*X';
                A2 = (A2.*exp( diag(1./(sum( A1,1) + eps) )*Nom*diag(1./(sum( X,2) + eps) ) ) ).^(1 + alphaSa);
                A2 = A2*diag(1./(sum(A2,1) + eps));
                
                Nom = (A1*A2)'*log(Y./(A1*A2*X + eps) + eps);
                X = (X.*exp(diag(1./(sum( A1*A2,1) + eps) )*Nom ) ).^(1 + alphaS);
                
                             
            end % type_alg
            
                
        else
               
          if type_alg == 1 % KL
                   Ap = A1;
               for t = 1:no_iter 
                   Nom = (Y./(A1*A2*X + eps))*(A2*X)' - alphaA1;
                   Nom(Nom <= 0) = eps;
                   A1 = (A1.*(Nom./( repmat((A2*X*ones(T,1))',m,1) + eps))).^(1 + alphaSa);
                   A1 = A1*diag(1./sum(A1,1));
                
               end
               for t = 1:no_iter 
                    Nom = A1'*(Y./(A1*A2*X + eps))*X' - alphaA2;
                    Nom(Nom <= 0) = eps;
                    A2 = (A2.*(Nom./(A1'*ones(m,1)*ones(1,T)*X' + eps))).^(1 + alphaSa);
                    A2 = A2*diag(1./sum(A2,1));
               end
               for t = 1:no_iter 
                   Nom = (A1*A2)'*(Y./(A1*A2*X + eps)) - alphaX;
                   Nom(Nom <= 0) = eps;
                   X = (X.*(Nom./(repmat((A1*A2)'*ones(m,1),1,T) + eps))).^(1 + alphaS);
               end
               
          elseif type_alg == 2 % Frobenius
               Ap = A1;
               for t = 1:no_iter 
                   Nom = Y*(A2*X)' - alphaA1;
                   Nom(Nom <= 0) = eps;
                   A1 = (A1.*(Nom./(A1*A2*X*(A2*X)' + eps))).^(1 + alphaSa);
                   A1 = A1*diag(1./sum(A1,1));
               end
               for t = 1:no_iter 
                     Nom = A1'*Y*X' - alphaA2;
                     Nom(Nom <= 0) = eps;
                     A2 = (A2.*(Nom./(A1'*A1*A2*X*X' + eps))).^(1 + alphaSa);
                     A2 = A2*diag(1./sum(A2,1));
               end
               for t = 1:no_iter 
                    Nom = (A1*A2)'*Y - alphaX;
                    Nom(Nom <= 0) = eps;
                    X = (X.*(Nom./((A1*A2)'*A1*A2*X + eps))).^(1 + alphaS);
               end
               
          elseif type_alg == 3 % Amari alpha-divergence
               Ap = A1;
               for t = 1:no_iter 
                    Nom = ((Y./(A1*A2*X + eps)).^alpha )*(A2*X)';
                    Nom(Nom <= 0) = eps;
                    A1 = (A1.*(Nom./( repmat((A2*X*ones(T,1))',m,1) + eps)).^(1/alpha) ).^(1 + alphaSa);
                    A1 = A1*diag(1./sum(A1,1));
               end
               for t = 1:no_iter 
                     Nom = A1'*((Y./(A1*A2*X + eps)).^alpha )*X';
                     Nom(Nom <= 0) = eps;
                     A2 = (A2.*(Nom./(A1'*ones(m,1)*ones(1,T)*X' + eps)).^(1/alpha) ).^(1 + alphaSa);
                     A2 = A2*diag(1./sum(A2,1));
               end
               for t = 1:no_iter 
                     Nom = (A1*A2)'*(Y./(A1*A2*X + eps)).^alpha;
                     Nom(Nom <= 0) = eps;
                     X = (X.*(Nom./(repmat((A1*A2)'*ones(m,1),1,T) + eps)).^(1/alpha) ).^(1 + alphaS);
               end
               
            elseif type_alg == 4 % SMART
               Ap = A1;
               for t = 1:no_iter 
                    Nom = ( log(Y./(A1*A2*X + eps) +eps) )*(A2*X)';
                    Nom(Nom <= 0) = eps;
                    A1 = (A1.*exp( Nom* diag(1./(sum( A2*X,2) +eps) ) ) ).^(1 + alphaSa);
                    A1 = A1*diag(1./(sum(A1,1) + eps));
               end
               for t = 1:no_iter 
                    Nom = A1'*( log(Y./(A1*A2*X + eps) + eps) )*X';
                    Nom(Nom <= 0) = eps;
                    A2 = (A2.*exp( diag(1./(sum( A1,1) + eps) )*Nom*diag(1./(sum( X,2) + eps) ) ) ).^(1 + alphaSa);
                    A2 = A2*diag(1./(sum(A2,1) + eps));
               end
               for t = 1:no_iter 
                    Nom = (A1*A2)'*log(Y./(A1*A2*X + eps) + eps);
                    Nom(Nom <= 0) = eps;
                    X = (X.*exp(diag(1./(sum( A1*A2,1) + eps) )*Nom ) ).^(1 + alphaS);
               end   
                             
          end % type_alg
                         
        end % if no_iter
 

        
                if (nr == max_restart)&(mod(k,50)==0)& (~mc_on | restart_mc_on)
                    norm_A = norm(abs(A1 - Ap),'fro');
                    fprintf(1, 'Restart %d,  %d-th alternating step\n',nr_best+1,k);
                end
        
                if sum(Index_norm)
                   if (nr == max_restart) & (((k < 50) & (mod(k,5)==0)) | ((k>49) & ((mod(k,50)==0)))) 
                       
                       s_dist = s_dist + 1;
                       k_select(s_dist) = k;
                       Z = A1*A2*X + eps;
                       Z = diag(1./(sqrt(var(Z')) + eps))*Z;
                      
                        dist_Fro(s_dist) = norm(Y - Z,'fro'); 
                        dist_KL(s_dist) = sum(sum(Y.*log(Y./Z + eps) - Y + Z)); 
                        dist_KL2(s_dist) = sum(sum(Z.*log(Z./Y + eps) + Y - Z)); 
                        dist_Pearson(s_dist) = sum(sum( ((Y - Z).^2)./Z ));
                        dist_Hellinger(s_dist) = sum(sum( (sqrt(Z) - sqrt(Y)).^2 )); 
                        dist_JS_rel(s_dist) = sum(sum(2*Y.*log(2*Y./(Y + Z) + eps) + Z - Y));  
                        dist_JS_rel2(s_dist) = sum(sum(2*Z.*log(2*Z./(Y + Z) + eps) - Z + Y));  
                        Zy = Y + Z; 
                        dist_JS(s_dist) = sum(sum(Y.*log(2*Y./Zy + eps) + Z.*log(2*Z./Zy + eps) ));  
                        dist_AG_rel(s_dist) = sum(sum(Zy.*log(.5*Zy./Y + eps) + Y - Z));  
                        dist_AG(s_dist) = sum(sum(.5*Zy.*log(.5*Zy./sqrt(Y.*Z) + eps)));  
                        dist_J(s_dist) = sum(sum( .5*(Y - Z).*log(Y./Z + eps) ));  
                        dist_Chi(s_dist) = sum(sum( ((Y + Z).*(Y - Z).^2)./(Y.*Z) ));  
                        dist_Tria(s_dist) = sum(sum( ((Y - Z).^2)./(Y + Z) ));  
                    end % if multiple
                end % if sum
                
            
end % while (k)
  
% Outer KL divergence
Z = AL*A1*A2*X; 
Z_outer = norm(Z,'fro') + eps;
KL_outer = sum(sum(Y_true.*log((Y_true + eps)./(Z + eps)) - Y_true + Z))/Z_outer;
         
          if (nr == 0) | (KL_outer < KL_outer_temp)
              A1_best = A1; A2_best = A2; X_best = X; KL_outer_temp = KL_outer; nr_best = nr;
           end % multi-conditions
           
   nr = nr + 1;
   
   if nr <=max_restart
      fprintf(1, 'Restart %d, Kullback-Leibler divergence = %e\n',	nr, KL_outer);
   end
   
end % while (restarts)

% One-Variance scaling
X(X <= 0) = eps;
A = A1*A2;

Distance_output = cell(length(s_dist),1);
Distance_output(1) = {[]};
Distance_output(2) = {[]};
Distance_output(3) = {[]};
Distance_output(4) = {[]};
Distance_output(5) = {[]};
Distance_output(6) = {[]};
Distance_output(7) = {[]};
Distance_output(8) = {[]};
Distance_output(9) = {[]};
Distance_output(10) = {[]};
Distance_output(11) = {[]};
Distance_output(12) = {[]};
Distance_output(13) = {[]};
Distance_output(14) = {[]};

if sum(Index_norm)
   Distance_output(1) = {k_select}; 
   Distance_output(2) = {dist_Fro};
   Distance_output(3) = {dist_KL};
   Distance_output(4) = {dist_KL2};
   Distance_output(5) = {dist_Pearson};
   Distance_output(6) = {dist_Hellinger};
   Distance_output(7) = {dist_JS_rel};
   Distance_output(8) = {dist_JS_rel2};
   Distance_output(9) = {dist_JS};
   Distance_output(10) = {dist_AG_rel};
   Distance_output(11) = {dist_AG};
   Distance_output(12) = {dist_J};
   Distance_output(13) = {dist_Chi};
   Distance_output(14) = {dist_Tria};
end

