% NMFLAB for Signal Processing written by A. Cichocki and R. Zdunek 
% in cooperation with other members of Laboratory for Advanced Brain Signal
% Processing, BSI, RIKEN, Saitama, JAPAN

function [A,X,Distance_output]=nmf_conv(Y,r,Index_norm, A_true, S, mc_on, NoAlts,restart_mc_on,AL,Y_true,type_alg_A, type_alg_X, max_restart, no_conv)
%
% Non-negative Matrix Factorization (NMF) with convolutive algorithms
%
% [A,X]=nmf_IP(Y, r, Index_norm, A_true, S, mc_on, NoAlts, restart_mc_on, AL, Y_true, type_alg_A, type_alg_X, ...
%       max_restart, no_iter, alpha0, tau) 
%       produces mixing matrix A of dimension [m by r],
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
%       > type_alg_A:    indicates the selected algorithm for computation
%                        of the mixing matrix,
%       > type_alg_X:    indicates the selected algorithm for computation of the sources,  
%       > max_restart:   number of restarts, 
%       > no_iter:       number of inner iterations, 
%       > alpha0:        initial magnitude in the exponential model for regularization parameter,
%       > tau:           damping factor in the exponential model for regularization parameter,
% 
% OUTPUTS:
%       > A:               estimated mixing matrix,
%       > X:               estimated source matrix,
%       > Distance_output: structures of different divergences measured between "Y" and estimated "AX" versus iterations,
%
% #########################################################################
A = [];
X = [];
if (nargin < 14) | isempty(no_conv) | (no_conv < 0) | max(size(no_conv) > 1)
   disp('Incorrect number of convolutive components');
   return
end
if (nargin < 13) | isempty(max_restart)  | (max_restart < 0) | max(size(max_restart) > 1)
   disp('Number of restarts must be given correctly');
   return
end
if (nargin < 12) | isempty(type_alg_X) | (type_alg_X < 1) | max(size(type_alg_X) > 1)
   disp('Incorrect algorithm for X');
   return
end
if (nargin < 11) | isempty(type_alg_A) | (type_alg_A < 1) | max(size(type_alg_A) > 1)
   disp('Incorrect algorithm for A');
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

[m,T]=size(Y);
niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)
niter_sample = 30; % maximum number of iterations for each random sample
epsil_normA = 1E-14; % tolerance for alternating

if mc_on & ~restart_mc_on
    max_restart = 0;
end
if ~isempty(NoAlts)
    niter_selected = NoAlts;
end

if (type_alg_A == 3) & (size(A_true,1) ~= size(Y,1))
   disp('Multilayer technique cannot be used with A fixed');
   Distance_output = [];
   return
end

% Declaration for A and X
A=zeros(m,r,no_conv+1);
Ap = A;
X=zeros(r,T);
Ainit = A;
Xinit = X;
Z = zeros(m,T);
KL_outer_temp = 0;
Z_outer = 0;
nr = 0; restart_on = 0; norm_A = 10; nr_best = -1;
m_sx = 1:m; r_sx = 1:r; T_sx = 1:T;
s_dist = 0;
delta = .1;

while (nr <= max_restart)
        
   % Initialize random A and X
       if ~nr & (~mc_on | restart_mc_on)   
          Ainit_1(m_sx',r_sx) = abs(repmat(.1*sin(2*pi*.1*m_sx'),1,r) + repmat(.1*cos(2*pi*.1*r_sx),m,1) + repmat(cos(2*pi*.471*m_sx'),1,r) + repmat(sin(2*pi*.471*r_sx),m,1));
          Ainit_1 = Ainit_1/max(max(Ainit_1));
          Ainit = repmat(Ainit_1,[1 1 no_conv+1]);
        
          Xinit(r_sx',T_sx) = abs(repmat(.1*sin(2*pi*.1*r_sx'),1,T) + repmat(.1*cos(2*pi*.1*T_sx),r,1) + repmat(cos(2*pi*.471*r_sx'),1,T) + repmat(sin(2*pi*.471*T_sx),r,1));
          Xinit = Xinit/max(max(Xinit));
       else
          Ainit=rand(m,r,no_conv+1);
          Xinit=rand(r,T);
       end
       
        % Normalization of initial guess
       Ainit = Ainit.*repmat((1./sum(Ainit,1)),[m,1,1]);
        
        if (nr == max_restart)&(max_restart > 0)
           A = A_best;
           X = X_best;
        else
           A = Ainit;
           X = Xinit;
        end % initial guess assignment
    
    Yx = zeros(m,T);
    n = 0; k = 0;

    
while ((k <= niter_sample)&(nr < max_restart)) | ((k <= niter_selected)&(nr == max_restart)&(norm_A > epsil_normA)& isempty(NoAlts)) | ((k <= niter_selected)&(nr == max_restart)& (NoAlts > 0)) 
 
k = k + 1;
           
l = mod(k-1,no_conv+1)+1;
%Xs = (circshift(X',l-1))';
Xs = X(:, (no_conv-l+2):(T-l+1));
Ys = Y(:, (no_conv-l+2):(T-l+1));

   
        switch type_alg_A
                            
           case 1 % ALS    
                          
              Ap = A; 
              Xpinv = Xs'*pinv(Xs*Xs');
              A(:,:,l) = max(1E6*eps, Ys*Xpinv);  
              A = A.*repmat((1./sum(A,1)),[m,1,1]);
                           
           case 2 % EMML    
              
              Ap = A; 
              A(:,:,l) = A(:,:,l).*((Y./(Yz + eps))*Xs')./repmat(sum(Xs,2)',m,1); 
              A = A.*repmat((1./sum(A,1)),[m,1,1]);    
              
           case 3 % Fixed A

              Ap = A; 
              niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)
              A = A_true + eps;    
              A = A.*repmat((1./sum(A,1)),[m,1,1]);

        end % type_alg_A 
   

        switch type_alg_X
                        
            case 1 % ALS
                
               Xs = max(1E6*eps,pinv(squeeze(A(:,:,l))'*squeeze(A(:,:,l)))*squeeze(A(:,:,l))'*Ys);
              
            case 2 % EMML
              
              Zy = (circshift((Y./(Yz + eps))',-l+1))';
               if l > 1
                  Zy(:,(T-l+2):T) = zeros(m,l-1);
               end  
              X = X.*(A(:,:,l)'*Zy)./repmat(sum(A(:,:,l),1)',1,T);  
                
            case 3 % Fixed X
            
              X = S + eps;   
              niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)

        end % type_alg_X
        
        
      %  X = (circshift(Xs',-l+1))';
      

                if (nr == max_restart)&(mod(k,50)==0)& (~mc_on | restart_mc_on)  
                   norm_A = sqrt(sum(sum(sum((A - Ap).^2))));
                   fprintf(1, 'Restart %d,  %d-th alternating step\n',nr_best+1,k);
                end
                                
                if sum(Index_norm)
                   if (nr == max_restart) & (((k < 50) & (mod(k,5)==0)) | ((k>49) & ((mod(k,50)==0)))) 
                       
                       s_dist = s_dist + 1;
                       k_select(s_dist) = k;
                       Z = A*X + eps;
                       Z = diag(1./sqrt(var(Z')))*Z;
                       Y = Y + eps;
                      
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
if ndims(A) > 2
   Yz = zeros(m,T);
   for t = 1:size(A,3)
       Ytmp = squeeze(A(:,:,t))*X;
       Yz = Yz + Ytmp;
   end
   Z = AL*Yz;
else
    Z = AL*A*X; 
end
Z_outer = norm(Z,'fro') + eps;
KL_outer = sum(sum(Y_true.*log((Y_true + eps)./(Z + eps)) - Y_true + Z))/Z_outer;
         
           if (nr == 0) | (KL_outer < KL_outer_temp)
              A_best = A; X_best = X; KL_outer_temp = KL_outer; nr_best = nr;
           end % multi-conditions
        
   nr = nr + 1;
   
   if nr <=max_restart
      fprintf(1, 'Restart %d, Kullback-Leibler divergence = %e\n',	nr, KL_outer);
   end
   
end % while (restarts)

% One-Variance scaling
X(X <= 0) = eps;

if ndims(A) > 2
   A =  sum(A,3);
end

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


 