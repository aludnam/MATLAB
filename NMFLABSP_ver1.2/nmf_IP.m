% NMFLAB for Signal Processing written by A. Cichocki and R. Zdunek 
% in cooperation with other members of Laboratory for Advanced Brain Signal
% Processing, BSI, RIKEN, Saitama, JAPAN

function [A,X,Distance_output]=nmf_IP(Y,r,Index_norm, A_true, S, mc_on, NoAlts,restart_mc_on,AL,Y_true,type_alg_A, ... 
type_alg_X, max_restart, no_iter, reg_oper_A, reg_oper_X, reg_param_A, reg_param_X, weightA, weightX, alphaA0, alphaX0, ...
tauA, tauX, alphaA_reg_const, alphaX_reg_const, alphaA, alphaX, p_norm, window_samples, overlapping_focuss, beta, niter_selected)
%
% Non-negative Matrix Factorization (NMF) with Interior-Point (IP) algorithms
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

if (nargin < 33) | isempty(niter_selected) | max(size(niter_selected) > 1)
   disp('Default number of alternating steps');
   niter_selected = 1000;
end
if (nargin < 32) | isempty(beta) | max(size(beta) > 1)
   disp('Incorrect beta parameter in NRALS');
   return
end
if (nargin < 31) | isempty(overlapping_focuss) | (overlapping_focuss < 0) | max(size(overlapping_focuss) > 1)
   disp('Incorrect overlapping');
   return
end
if (nargin < 30) | isempty(window_samples) | (window_samples < 1) | max(size(window_samples) > 1)
   disp('Incorrect number of samples in the window');
   return
end
if (nargin < 29) | isempty(p_norm) | max(size(p_norm) > 1)
   disp('Incorrect parameter p in M-FOCUSS');
   return
end
if (nargin < 28) | isempty(alphaX) | max(size(alphaX) > 1)
   disp('Incorrect the sparsity parameter for X');
   return
end
if (nargin < 27) | isempty(alphaA) | max(size(alphaA) > 1)
   disp('Incorrect the sparsity parameter for A');
   return
end
if (nargin < 26) | isempty(alphaX_reg_const) | (alphaX_reg_const < 0) | max(size(alphaX_reg_const) > 1)
   disp('Incorrect the constant regularization parameter for X');
   return
end
if (nargin < 25) | isempty(alphaA_reg_const) | (alphaA_reg_const < 0) | max(size(alphaA_reg_const) > 1)
   disp('Incorrect the constant regularization parameter for A');
   return
end
if (nargin < 24) | isempty(tauX) | max(size(tauX) > 1)
   disp('Incorrect the damping factor in the exponential model for X');
   return
end
if (nargin < 23) | isempty(tauA) | max(size(tauA) > 1)
   disp('Incorrect the damping factor in the exponential model for A');
   return
end
if (nargin < 22) | isempty(alphaX0) | max(size(alphaX0) > 1)
   disp('Incorrect the initial magnitude in the exponential model for X');
   return
end
if (nargin < 21) | isempty(alphaA0) | max(size(alphaA0) > 1)
   disp('Incorrect the initial magnitude in the exponential model for A');
   return
end
if (nargin < 20) | isempty(weightX) | (weightX < 1) | max(size(weightX) > 1)
   disp('Incorrect weighting in X');
   return
end
if (nargin < 19) | isempty(weightA) | (weightA < 1) | max(size(weightA) > 1)
   disp('Incorrect weighting in A');
   return
end
if (nargin < 18) | isempty(reg_param_X) | (reg_param_X < 1) | max(size(reg_param_X) > 1)
   disp('Incorrect type of regularization parameter in X');
   return
end
if (nargin < 17) | isempty(reg_param_A) | (reg_param_A < 1) | max(size(reg_param_A) > 1)
   disp('Incorrect type of regularization parameter in A');
   return
end
if (nargin < 16) | isempty(reg_oper_X) | (reg_oper_X < 1) | max(size(reg_oper_X) > 1)
   disp('Incorrect type of regularization operator in X');
   return
end
if (nargin < 15) | isempty(reg_oper_A) | (reg_oper_A < 1) | max(size(reg_oper_A) > 1)
   disp('Incorrect type of regularization operator in A');
   return
end
if (nargin < 14) | isempty(no_iter) | (no_iter < 1) | max(size(no_iter) > 1)
   disp('Incorrect number of inner iterations');
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
% if min(sum(Y,2)) == 0
%     disp('Not all entries in a row can be zero');
%     return
% end

[m,T]=size(Y);
%niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)
niter_sample = 30; % maximum number of iterations for each random sample
epsil_normA = 1E-5; % tolerance for alternating

if mc_on & ~restart_mc_on
    max_restart = 0;
end
if ~isempty(NoAlts)
    niter_selected = NoAlts;
end

if (type_alg_A == 7) & (size(A_true,1) ~= size(Y,1))
   disp('Multilayer technique cannot be used with A fixed');
   Distance_output = [];
   return
end

% M-FOCUSS
n_overlap = overlapping_focuss*window_samples/100;
if overlapping_focuss
    t_max = ceil((T - n_overlap)/(window_samples - n_overlap));
else
    t_max = floor(T/window_samples);
end

% Declaration for A and X
A=zeros(m,r);
Ap = A;
X=zeros(r,T);
Ainit = A;
Xinit = X;
Z = zeros(m,T);
E = ones(r);
KL_outer_temp = 0;
Z_outer = 0;
nr = 0; restart_on = 0; norm_A = 10; nr_best = -1;
m_sx = 1:m; r_sx = 1:r; T_sx = 1:T;
s_dist = 0;
delta = .1;

while (nr <= max_restart)
        
   % Initialize random A and X
       if ~nr & (~mc_on | restart_mc_on)   
          Ainit(m_sx',r_sx) = abs(repmat(.1*sin(2*pi*.1*m_sx'),1,r) + repmat(.1*cos(2*pi*.1*r_sx),m,1) + repmat(cos(2*pi*.471*m_sx'),1,r) + repmat(sin(2*pi*.471*r_sx),m,1));
          Ainit = Ainit/max(max(Ainit));
        
          Xinit(r_sx',T_sx) = abs(repmat(.1*sin(2*pi*.1*r_sx'),1,T) + repmat(.1*cos(2*pi*.1*T_sx),r,1) + repmat(cos(2*pi*.471*r_sx'),1,T) + repmat(sin(2*pi*.471*T_sx),r,1));
          Xinit = Xinit/max(max(Xinit));
       else
          Ainit=rand(m,r);
          Xinit=rand(r,T);
       end
       
        % Normalization of initial guess
       Ainit = Ainit*diag(1./sum(Ainit,1));
        
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

switch reg_oper_X

    case 1 % Ones
        
         CX = ones(r); 
         CF = ones(m); 
         
    case 2 % Identity
        
         CX = eye(r); 
         CF = eye(m); 
                
    case 3 % Diff L1
        
         CX = get_l(r,1);
         CX = CX'*CX;
         CF = get_l(m,1);
         CF = CF'*CF;
                  
    case 4 % Diff L2
        
         CX = get_l(r,2);
         CX = CX'*CX;
         CF = get_l(m,2);
         CF = CF'*CF;     
         
    case 5 % PCA-Cov
        
        Cr = cov((Y - A*X)');
        [V,D] = eig(Cr);
        V = fliplr(V); D = flipud(fliplr(D));
        Vr = V(:,1:r); Er = Vr'*Cr;
        CX = Er*Er'/T;
        CF = V'*Cr*Cr'*V/T; 
end


switch reg_param_X

    case 1 % Zero
        
         alphaX_reg = 0;
         
    case 2 % Exp
        
         alphaX_reg = alphaX0*exp(-k/tauX);
        
    case 3 % Const
        
         alphaX_reg = alphaX_reg_const;
         
    case 4 % Alpha_reg res(k)
        
         e_res = 0.5*norm(Y - A*X,'fro');
         alphaX_reg = alphaX_reg_const*e_res;
                
    case 5 % EM-PCA
        
        sigma_v = trace((Y - A*X)*Y')/(m*T);
        alphaX_reg = sigma_v;
end

switch weightX
    
    case 1 % Identity
        
        W = eye(m);
        
    case 2 % BLUE
        
        W = pinv(cov((Y - A*X)'));  
end


        for t = 1:no_iter % inner iterations

                  switch type_alg_X

                      case 1 % ALS

                         X = max(1E6*eps,pinv(A'*A)*A'*Y);    
                                               
                      case 2 % Sparse ALS

                         X = max(1E6*eps,pinv(A'*A)*(A'*Y - alphaX)); 
                                                  
                      case 3 % BLUE ALS
                     
                         W = pinv(cov((Y - A*X)')); 
                         X = max(1E6*eps,pinv(A'*W*A)*A'*W*Y);       

                      case 4 % Regularized weighted sparse ALS
                                         
                         X = max(1E6*eps,pinv((A'*W*A) +  alphaX_reg*CX)*(A'*W*Y - alphaX));      
                                                  
                      case 5 % Self-compensated regularized sparse ALS
                          
                         X = max(1E6*eps,pinv(A'*A +  alphaX_reg*CX)*(A'*Y - alphaX + alphaX_reg*CX*X));   
                         
                      case 6 % ROM-FOCUSS
                 
                            s_z = [];
                            for t = 1:t_max

                                if overlapping_focuss & (t < t_max)
                                   s_z = floor((t-1)*(window_samples - n_overlap))+1:floor(t*window_samples - (t-1)*n_overlap); 
                                elseif  overlapping_focuss & (t == t_max)
                                   s_z = floor((t-1)*(window_samples - n_overlap))+1:T;  
                                else
                                   s_z = window_samples*(t-1)+1:window_samples*t;
                                end

                                cx = sqrt(sum(X(:,s_z).^2,2));
                                Dx = diag(abs(cx).^(2 - p_norm));
                                X(:,s_z) = max(1E6*eps,Dx*A'*pinv(A*Dx*A' + alphaX_reg*CF)*Y(:,s_z)); 

                            end % t

                      case 7 % MALS by Wang-Hopke
                          
                         CX = min(diag(sum(abs(Y'*A))./(sum(X')+eps)),diag(sum(A'*A))); 
                         X = max(1E6*eps,pinv(A'*A +  alphaX_reg*CX)*(A'*Y - alphaX + alphaX_reg*CX*X));    
                                                  
                      case 8 % Hierarchical ALS
                                            
                         for j = 1:r 
                             Ys = Y - A*X + A(:,j)*X(j,:);
                             X(j,:) = max(1E6*eps,(A(:,j)'*Ys - alphaX*sign(X(j,:)))/(A(:,j)'*A(:,j)));      
                         end
                                           
                         
                      case 9 % Nonlinear RALS
                          
                          sigma = norm(Y - A*X,'fro')/(r*T);
                          X = max(1E6*eps,inv(A'*A + sigma*eye(r))*(A'*Y - sigma*tanh(beta*X)));
                         
                      case 10 % Fixed X

                         X = S + eps;   
                         niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)

                  end % type_alg_X

        end % t for X

        
  switch reg_oper_A

    case 1 % Ones
        
         CA = ones(r); 
         
    case 2 % Identity
        
         CA = eye(r); 
        
    case 3 % Diff L1
        
         CA = get_l(r,1);
         CA = CA'*CA;
         
    case 4 % Diff L2
        
         CA = get_l(r,2);
         CA = CA'*CA;
                
    case 5 % PCA-Cov
        
        Cr = cov((Y - A*X)');
        [V,D] = eig(Cr);
        V = fliplr(V); D = flipud(fliplr(D));
        Vr = V(:,1:r); Er = Vr'*Cr;
        CA = Er*Er'/T; 
  end
      

 switch reg_param_X

    case 1 % Zero
        
         alphaA_reg = 0;
         
    case 2 % Exp
        
         alphaA_reg = alphaA0*exp(-k/tauA);
        
    case 3 % Const
        
         alphaA_reg = alphaA_reg_const;
         
    case 4 % Alpha_reg res(k)
        
         e_res = 0.5*norm(Y - A*X,'fro');
         alphaA_reg = alphaA_reg_const*e_res;
                
    case 5 % EM-PCA
        
        sigma_v = trace((Y - A*X)*Y')/(m*T);
        alphaA_reg = m*sigma_v;
 end
  
 
switch weightA
    
    case 1 % Identity
        
        W = eye(m);
        
    case 2 % BLUE
        
        W = cov((Y - A*X)');  
end

        
        for t = 1:no_iter % inner iterations

                  switch type_alg_A

                      case 1 % ALS

                         Ap = A;        
                         A = max(1E6*eps, Y*X'*pinv(X*X'));  
                         A = A*diag(1./sum(A,1));   

                      case 2 % Sparse ALS   
                          
                         Ap = A;        
                         A = max(1E6*eps, (Y*X' - alphaA)*pinv(X*X'));  
                         A = A*diag(1./sum(A,1));    
                                          
                      case 3 % Regularized weighted sparse ALS 

                         Ap = A;        
                         A = max(1E6*eps, (Y*X' - alphaA*W*ones(m,r))*pinv((X*X') + alphaA_reg*CA));  
                         A = A*diag(1./sum(A,1));      
                      
                      case 4 % Self-compensated regularized sparse ALS

                         Ap = A;        
                         A = max(1E6*eps, (Y*X' - alphaA + alphaA_reg*A*CA)*pinv(X*X' + alphaA_reg*CA));  
                         A = A*diag(1./sum(A,1));      
                                                     
                      case 5 % MALS by Wang-Hopke   
                         
                         Ap = A;
                         CA = min(diag(sum(abs(Y*X'))./(sum(A)+eps)),diag(sum(X*X')));
                         A = max(1E6*eps, (Y*X' - alphaA + alphaA_reg*A*CA)*pinv(X*X' + alphaA_reg*CA));  
                         A = A*diag(1./sum(A,1));   
                      
                      case 6 % Hierarchical ALS
                         
                         Ap = A;
                         for j = 1:r
                             Ys = Y - A*X + A(:,j)*X(j,:);
                             A(:,j) = max(1E6*eps, Ys*X(j,:)'/(X(j,:)*X(j,:)'));  
                         end
                         A = A*diag(1./sum(A,1));      
                         
                      case 7 % Fixed A

                         niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)
                         A = A_true + eps;    
                         A = A*diag(1./(sum(A,1) + eps));

                  end % type_alg_A 
                              
        end % t for A
                           
         
                if (nr == max_restart)&(mod(k,50)==0)& (~mc_on | restart_mc_on)  
                   norm_A = norm(abs(A - Ap),'fro');
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

Z = AL*A*X; 
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
%X = max(1E6*eps,pinv(A'*A)*A'*Y); 

%X(X <= 0) = eps;

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



 