function [u,covar,t,iter] = fit_mix_2D_gaussian( X,M )
%
% fit_mix_2D_gaussian - fit parameters for a 2D mixed-gaussian distribution using EM algorithm
%
% format:   [u,covar,t,iter] = fit_mix_2D_gaussian( X,M )
%
% input:    X   - input samples, Nx2 vector
%           M   - number of gaussians which are assumed to compose the distribution
%
% output:   u       - fitted mean for each gaussian (each mean is a 2x1 vector)
%           covar   - fitted covariance for each gaussian. this is a 2x2xM matrix.
%           t       - probability of each gaussian in the complete distribution
%           iter    - number of iterations done by the function
%

% run with default values
if ~nargin
    M   = 1 + floor(rand*5);
    for m = 1:M
        D               = diag( rand(1,2)*10 );         % create a leagal cov matrix
        A               = randn(2,2);                   % ...first choose the eigenvalues
        A               = A/det(A);                     % ...and later, a random transform
        covar(:,:,m)    = A*D*A';                       % ...the result -> cov matrix
        u(:,m)          = randn(2,1)*5;                 % the mean of the gaussians
        prob(m)         = rand;                         % each gaussian probability
    end
    X = build_mix_2D_gaussian( u,covar,prob,5000 );     % create the distribution samples
    [u,covar,t,iter] = fit_mix_2D_gaussian( X,M );      % estimate it...
    return
end

% set X dimensions
if (size(X,2)~=2),
    X = X.';
end

% initialize and initial guesses
N           = size( X,1 );
Z           = ones(N,M) * 1/M;                      % indicators vector
P           = zeros(N,M);                           % probabilities vector for each sample and each model
t           = ones(1,M) * 1/M;                      % distribution of the gaussian models in the samples
u           = [linspace(min(X(:,1)),max(X(:,1)),M);...
        linspace(min(X(:,2)),max(X(:,2)),M)];       % mean vector
covar       = repmat( cov(X),[1,1,M] )/ sqrt(M);    %  vector of covariance matrices
C           = 1/2/pi;                               % just a constant
Ic          = ones(N,1);                            % - enable a row replication by the * operator
Ir          = ones(1,2);                            % - enable a column replication by the * operator
Q           = zeros(N,M);                           % user variable to determine when we have converged to a steady solution
thresh      = 1e-3;
step        = N;
last_step   = inf;
iter        = 0;
min_iter    = 10;

% main convergence loop, assume gaussians are 1D
while ((( abs((step/last_step)-1) > thresh) & (step>(N*eps)) ) | (iter<min_iter) ) 
    
    % E step
    % ========
    % estimate the probability matrix:
    % P(n,m) = C / d * exp( -1/2 * ( (X(n,:)-U.')*invS*(X(n,:)'-U) ) );
    % calc each column at once, therefore, inv(S) = [s22 -s12;-s21 s11]/det(S)
    % which will enable a vectoric X of dimension 2 to be calculated in
    % (X-U)'*inv(S)*(X-U) => 1D vector
    Q = Z;
    for m = 1:M
        S       = covar(:,:,m);
        d       = det(S);
        U       = u(:,m);
        P(:,m)  = C ./ (Ic*d) .* exp( -1/2/d * (...
            S(2,2)*(X(:,1)-U(1)).^2 + S(1,1)*(X(:,2)-U(2)).^2 - ...
            (S(2,1)+S(1,2))*(X(:,1)-U(1)).*(X(:,2)-U(2)) ) );
        Z(:,m)  = (P(:,m)*t(m))./(P*t(:));
    end
        
    % estimate convergence step size and update iteration number
    prog_text   = sprintf(repmat( '\b',1,(iter>0)*12+ceil(log10(iter+1)) ));
    iter        = iter + 1;
    last_step   = step * (1 + eps) + eps;
    step        = sum(sum(abs(Q-Z)));
    fprintf( '%s%d iterations\n',prog_text,iter );

    % M step
    % ========
    Zm              = sum(Z);               % sum each column
    Zm(find(Zm==0)) = eps;                  % avoid devision by zero
    t               = Zm/N;
    for m = 1:M
        dist        = X - Ic*(u(:,m)');
        covar(:,:,m)= dist' * ( (Z(:,m)*Ir) .* dist )/Zm(m);    % covariance estimation
        u(:,m)      = (X')*Z(:,m) / Zm(m);                      % mean estimation
    end
end

% plot the fitted distribution
% =============================
plot_mix_gaussian( u,covar,t );
