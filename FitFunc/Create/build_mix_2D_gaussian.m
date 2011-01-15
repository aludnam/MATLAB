function [X,N] = build_mix_2D_gaussian( u_c,covar_c,p_c,N )
%
% build_mix_2D_gaussian - build a distribution of a mixed gaussian 
% 
% format:   [X,N] = build_mix_2D_gaussian( u_c,covar_c,p_c,N )
%
% input:    u_c     - 2xM matrix, mean of the gaussians in the distribution
%           covar_c - 2x2xM matirx, covariance of the gaussians
%           p_c     - 1xM vector, ratio between the gaussians in the distribution
%           N       - total number of samples, this number might change a little
%
% output:   X       - the random variable matrix Nx2
%           N       - the length of the given 2D vector
%

% initialize
M           = length(p_c);                  % the number of mixed gaussians
p_c         = round( p_c / sum(p_c) * N );  % number of samples for each gaussian
N           = sum( p_c );                   % total number of samples
X           = zeros(N,2);                   % initialize output samples 2D vector
idx         = 1;                            % counter
actual_cov  = zeros(2,2,M);                 % the actual covariance calculated from the samples
actual_mean = zeros(2,M);                   % the actual mean calculated from the samples

% build distribution
for m = 1:M
    if (idx>1)
        span = sum(p_c(1:idx-1)):sum(p_c(1:idx));
    else
        span = 1:p_c(1);
    end
    X(span,:) = mvnrnd( u_c(:,idx),covar_c(:,:,idx),length( span ) );
    idx       = idx + 1;
    actual_cov(:,:,m) = cov( X(span,:) );
    actual_mean(:,m)  = mean( X(span,:) )';
end

% plot the ACTUAL distribution
plot_mix_gaussian( actual_mean,actual_cov,p_c/sum(p_c),X );
