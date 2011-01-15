function result = fit_ML_normal( x,hAx )
% fit_ML_normal - Maximum Likelihood fit of the normal distribution of i.i.d. samples!.
%                  Given the samples of a normal distribution, the PDF parameter is found
%
%    fits data to the probability of the form: 
%        p(r) = sqrt(1/2/pi/sig^2)*exp(-((r-u)^2)/(2*sig^2))
%    with parameters: u,sig^2
%
% format:   result = fit_ML_normal( x,hAx )
%
% input:    x   - vector, samples with normal distribution to be parameterized
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      sig^2,u          - fitted parameters
%                      CRB_sig2,CRB_u   - Cram?r-Rao Bound for the estimator value
%                      RMS              - RMS error of the estimation 
%                      type             - 'ML'
%
% example:  fit_ML_normal( randn(1,10000)*3 - 1 )
%               or
%           figure; 
%           fit_ML_normal( randn(1,10000)*3 - 1,gca );
%

%
% Algorithm
% ===========
%
% We use the ML algorithm to estimate the PDF from the samples.
% The normal destribution is given by:
%
%    p(x;u,sig^2) = sqrt(1/2/pi/sig^2)*exp(-((x-u).^2)/(2*sig^2))
%
%    where x are the samples which distribute by the function p(x;u,sig^2)
%            and are assumed to be i.i.d !!!
%
% The ML estimator is given by:
%
%    a         = parameters vector = [u,sig^2]
%    f(Xn,a)   = sqrt(1/2/pi/sig^2)*exp(-((Xn-u).^2)/(2*sig^2))
%    L(a)      = f(X,a) = product_by_n( f(Xn,a) )
%              = (2*pi*sig^2)^(-N/2)*exp(-sum((Xn-u)^2)/(2*sig^2))
%    log(L(a)) = -N/2*log(2*pi) - N/2*log(sig^2) - sum((Xn-u)^2)/(2*sig^2)
%
%    The maximum likelihood point is found by the derivative of log(L(a)) with respect to "a":
%
%    diff(log(L(a)),u) = -2*sum(Xn-u)/(2*sig^2) = N/sig^2 * ( u - sum(Xn)/N )
%                    = J(u) * (u_estimation - u)  
%    diff(log(L(a)),sig^2) = -N/(2*sig^2) + sum((Xn-u)^2)/(2*sig^4)
%                          = N/(2*sig^4) * ( sum((Xn-u)^2)/N - sig^2 )
%                          = J(sig^2) * (sig^2_estimator - sig^2)
%
%    Therefore, the (efficient) estimators are given by:
%
%               u     = sum(Xn)/N
%               sig^2 = sum((Xn-u)^2)/N
%
%    The Cram?r-Rao Bounds for these estimator are:
%
%               VAR( u )     = 1/J(u)     = (sig^2) / N
%               VAR( sig^2 ) = 1/J(sig^2) = (2*sig^4) / N
%
%    NOTE: the ML estimator does not detect a deviation from the model.
%          therefore, check the RMS value !
%

if (nargin<1)
    error( 'fit_ML_normal - insufficient input arguments' );
end

% Estimation
% =============
x       = x(:);                 % should be column vectors !
N       = length(x);
u       = sum(x)/N;
sig2    = (x-u)'*(x-u)/N;
CRB_u   = sig2 / N;
CRB_sig2= (2*sig2^2) / N;
[n,x_c] = hist( x,100 );
n       = n / sum(n*abs(x_c(2)-x_c(1)));
y       = sqrt(1/2/pi/sig2)*exp(-((x_c-u).^2)/(2*sig2));
RMS     = sqrt( (y-n)*((y-n)')/ (x_c(2)-x_c(1))^2 / (length(x_c)-1) );

% finish summarizing results
% ============================
result = struct( 'u',u,'sig2',sig2,'CRB_u',CRB_u,'CRB_sig2',CRB_sig2,'RMS',RMS,'type','ML' );

% plot distribution if asked for
% ===============================
if (nargin>1)
    xspan = linspace(min(x),max(x),100);
    if ishandle( hAx )
        plot_normal( xspan,result,hAx,1 );
    else
        figure;
        plot_normal( xspan,result,gca,1 );
    end
end
