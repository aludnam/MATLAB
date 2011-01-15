function [A,X,Distance_output]=nmf_alpha_alg(Y,r,Index_norm,A_true,S, mc_on, NoAlts,restart_mc_on,AL,Y_true,type_alg_A,type_alg_X, max_restart,no_iter, alphaX, alphaS, beta_dec,beta,loss_fun)
%
% Non-negative Matrix Factorization (NMF) with multi-algorithms
%
% [A,X]=nmf_multi(Y,r,n_restart,gamma) produces matrix A of dimension [m by r]
% and matrix X of dimension [r by T] which are respectively an estimate of the mixture matrix and 
% an estimate of the source signals of the linear model: AX = Y, where Y is
% an observation matrix [m by T]. 
% Note: > m: sensor number.
%       > r: source number.
%       > T: sample number.
%       > max_restart: number of restarts in selection of an initial guess for A and X (by default 10)
% #########################################################################
A = [];
X = [];
if (nargin < 19) | isempty(loss_fun) | (loss_fun < 0) | max(size(loss_fun) > 1)
   disp('Incorrect the loss function in the projected gradient algorithm');
   return
end
if (nargin < 18) | isempty(beta) | max(size(beta) > 1)
   disp('Incorrect "alpha" in Amari alpha-divergence');
   return
end
if (nargin < 17) | isempty(beta_dec) | (beta_dec < 0) | max(size(beta_dec) > 1)
   disp('Incorrect "beta" in the Armijo rule');
   return
end
if (nargin < 16) | isempty(alphaS) | max(size(alphaS) > 1)
   disp('Incorrect parameter alphaS');
   return
end
if (nargin < 15) | isempty(alphaX) | max(size(alphaX) > 1)
   disp('Incorrect regularization parameter');
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
if min(sum(Y,2)) == 0
    disp('Not all entries in a row can be zero');
    return
end

Y = Y + eps;

[m,T]=size(Y);
niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)
niter_sample = 30; % maximum number of iterations for each random sample
epsil_normA = 1E-4; % tolerance for alternating
tol = 0.001; % Tolerance in the PG algorithm
gamma_NG = .1;

% Monte Carlo and alternatings adjustment
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

% Declaration for A and X
A=zeros(m,r);
Ap = A;
X=zeros(r,T);
Ainit = A;
Xinit = X;
nr_best = -1;
Z = zeros(m,T);
KL_outer_temp = 0;
nr = 0; restart_on = 0; norm_A = 10;
m_sx = 1:m; r_sx = 1:r; T_sx = 1:T; s_dist = 0;


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
 
    if (type_alg_X == 1) | (type_alg_A == 1)
       gradA = A*(X*X') - Y*X'; gradX = (A'*A)*X - A'*Y;
       initgrad = norm([gradA; gradX'],'fro');
       tolA = max(0.001,tol)*initgrad; tolX = tolA;
    end
   
    
while ((k <= niter_sample)&(nr < max_restart)) | ((k <= niter_selected)&(nr == max_restart)&(norm_A > epsil_normA)& isempty(NoAlts)) | ((k <= niter_selected)&(nr == max_restart)& (NoAlts > 0)) 
 
k = k + 1;
 
    switch type_alg_X
                           
                      
        case 1 % Projected gradient method
                      
            [X,gradX,iterX] = nmf_proj_grad(Y,A,X,tolX,no_iter,beta_dec,loss_fun,beta);
             if iterX==1,
                tolX = 0.1 * tolX; 
             end
        
        case 2 % IPG
                      
             X = ipg(A, Y, X, no_iter);
             
        case 3 % MRNSD
            
            if alphaX
               X = mrnsd_reg(A,Y,X,no_iter,[],[],0,alphaX);
            else
               X = mrnsd(A,Y,X,no_iter,[],[],0); 
            end    
                       
        case 4 % Relaxed Hellinger
            
            for t = 1:no_iter
                X = (X.*(.25 + 0.25*sqrt(1 + 8*A'*(Y./(A*X + eps))))).^(1+alphaS);
            end
        
        case 5 % Csiszar
            
            for t = 1:no_iter
                X = (X - A'*(A*X - Y)./repmat(r*diag(A'*A),1,T)).^(1+alphaS);
                X(X<=0) = eps;
            end
               
                     
        case 6 % Pinv
            
             X = max(1E6*eps,pinv(A)*Y);   
                         
        case 7 % Fixed X
            
            X = S + eps; 
            niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)    
                 
    end % switch for X
    
 switch type_alg_A
           
                  
        case 1 % Projected gradient method
            
            Ap = A;    
            [A,gradA,iterA] = nmf_proj_grad(Y',X',A',tolA,no_iter,beta_dec,loss_fun,beta); A = A'; gradA = gradA';
            if iterA==1,
               tolA = 0.1 * tolA;
            end

        case 2 % IPG
            
            Ap = A;    
            A = ipg(X', Y', A', no_iter); A = A'; 
       
        case 3 % MRNSD
            
            Ap = A;    
            A = mrnsd(X',Y',A',no_iter,[],[],0);
            A = A';
            
        case 4 % Relaxed Hellinger
            
            Ap = A;    
            for t = 1:no_iter
                A = A.*(.25 + .25*sqrt(1 + 8*(Y./(A*X + eps))*X'));
            end
            A = A*diag(1./(sum(A,1) + eps));
        
        case 5 % Csiszar
            
            Ap = A;    
            for t = 1:no_iter
                A = A - ((A*X - Y)*X')./repmat(m*diag(X*X')',m,1);
                A(A<=0) = eps;
            end
               
        case 6 % Pinv
            
                 Ap = A;        
                 A = max(1E6*eps,Y*pinv(X')');  
                 A = A*diag(1./sum(A,1));      
                           
        case 7 % Fixed A
            
            niter_selected = 1000;     % maximum number of iterations for the selected sample (can be adjusted)
            A = A_true + eps;    
            
    end % switch for A
     


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

X(X <= 0) = eps;

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


