function result = fit_ML_laplace( x,hAx )
% fit_ML_normal - Maximum Likelihood fit of the laplace distribution of i.i.d. samples!.
%                  Given the samples of a laplace distribution, the PDF parameter is found
%
%    fits data to the probability of the form: 
%        p(x) = 1/(2*b)*exp(-abs(x-u)/b)
%    with parameters: u,b
%
% format:   result = fit_ML_laplace( x,hAx )
%
% input:    x   - vector, samples with laplace distribution to be parameterized
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      u,b    - fitted parameters
%                      CRB_b  - Cram?r-Rao Bound for the estimator value
%                      RMS    - RMS error of the estimation 
%                      type   - 'ML'
%

%
% Algorithm
% ===========
%
% We use the ML algorithm to estimate the PDF from the samples.
% The laplace destribution is given by:
%
%    p(x;u,b) = 1/(2*b)*exp(-abs(x-u)/b)
%
%    where x are the samples which distribute by the function p(x;u,b)
%            and are assumed to be i.i.d !!!
%
% The ML estimator is given by:
%
%    a         = parameters vector = [u,b]
%    f(Xn,a)   = 1/(2*b)*exp(-abs(Xn-u)/b)
%    L(a)      = f(X,a) = product_by_n( f(Xn,a) )
%              = (2*b)^(-N) * exp( - sum( abs(Xn-u) )/b )
%    log(L(a)) = -N*log(2*b) - sum( abs(Xn-u) )/b
%
%    The maximum likelihood point is found by the derivative of log(L(a)) with respect to "a":
%
%    diff(log(L(a)),b)  = N/(b^2) * ( sum( abs(Xn-u) )/N - b )
%                       = J(b) * (b_estimation - b)  
%    diff(log(L(a)),m)  = (1/b) * sum( diff( abs(Xn-u),u ) ) => can't obtain a derivative
%                       But, u is the mean of the distribution, and therefore => u = mean(Xn)
%                       
%
%    Therefore, the (efficient) estimators are given by:
%
%               u = sum( Xn )/N
%               b = sum( abs(Xn-u) )/N
%
%    The Cram?r-Rao Bounds for these estimator are:
%
%               VAR( m ) = ?
%               VAR( b ) = 1/J(b) = b^2 / N
%
%    NOTE: the ML estimator does not detect a deviation from the model.
%          therefore, check the RMS value !
%

if (nargin<1)
    error( 'fit_ML_laplace - insufficient input arguments' );
end

% Estimation
% =============
x       = x(:);                 % should be column vectors !
N       = length(x);
u       = sum( x )/N;
b       = sum(abs(x-u))/N;
CRB_b   = b^2 / N;
[n,x_c] = hist( x,100 );
n       = n / sum(n*abs(x_c(2)-x_c(1)));
y       = 1/(2*b)*exp(-abs(x_c-u)/b);
RMS     = sqrt( (y-n)*((y-n)')/ (x_c(2)-x_c(1))^2 / (length(x_c)-1) );

% finish summarizing results
% ============================
result = struct( 'u',u,'b',b,'CRB_b',CRB_b,'RMS',RMS,'type','ML' );

% plot distribution if asked for
% ===============================
if (nargin>1)
    xspan = linspace(min(x),max(x),100);
    if ishandle( hAx )
        plot_laplace( xspan,result,hAx,1 );
    else
        figure;
        plot_laplace( xspan,result,gca,1 );
    end
end
