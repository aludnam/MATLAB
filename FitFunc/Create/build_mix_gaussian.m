function [X,N] = build_mix_gaussian( u_c,sig_c,p_c,N )
%
% build_mix_gaussian - build a distribution of a mixed gaussian 
% 
% format:   [X,N] = build_mix_gaussian( u_c,sig_c,p_c,N )
%
% input:    u_c     - vector, mean of the gaussians in the distribution
%           sig_c   - vector, Standard Deviation of the gaussians
%           p_c     - vector, ratio between the gaussians in the distribution
%           N       - total number of samples, this number might change a little
%
% output:   X       - the random variable vector
%           N       - the length of the given vector
%

% initialize
M   = length(p_c);      % the number of mixed gaussians
p_c = round( p_c / sum(p_c) * N );
N   = sum( p_c );
X   = zeros(N,1);
idx = 1;

% build distribution
for m = 1:M
    if (idx>1)
        span = sum(p_c(1:idx-1)):sum(p_c(1:idx));
    else
        span = 1:p_c(1);
    end
    X( span ) = randn( length(span),1 )*sig_c(idx) + u_c(idx);
    idx       = idx + 1;
end

% plot the distribution
plot_mix_gaussian( u_c,sig_c,p_c/sum(p_c),X );