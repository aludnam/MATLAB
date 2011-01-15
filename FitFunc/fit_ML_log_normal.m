function result = fit_ML_log_normal( x,hAx )
% fit_ML_normal - Maximum Likelihood fit of the log-normal distribution of i.i.d. samples!.
%                  Given the samples of a log-normal distribution, the PDF parameter is found
%
%    fits data to the probability of the form: 
%        p(x) = sqrt(1/(2*pi))/(s*x)*exp(- (log(x-m)^2)/(2*s^2))
%    with parameters: m,s
%
% format:   result = fit_ML_log_normal( x,hAx )
%
% input:    x   - vector, samples with log-normal distribution to be parameterized
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      m,s          - fitted parameters
%                      CRB_m,CRB_s  - Cram?r-Rao Bound for the estimator value
%                      RMS          - RMS error of the estimation 
%                      type         - 'ML'
%

%
% Algorithm
% ===========
%
% We use the ML algorithm to estimate the PDF from the samples.
% The log-normal destribution is given by:
%
%    p(x;m,s) = sqrt(1/(2*pi))/(s*x)*exp(- (log(x-m)^2)/(2*s^2))
%
%    where x are the samples which distribute by the function p(x;m,s)
%            and are assumed to be i.i.d !!!
%
% The ML estimator is given by:
%
%    a         = parameters vector = [m,s]
%    f(Xn,a)   = sqrt(1/(2*pi))/(s*Xn)*exp(- (log(Xn-m)^2)/(2*s^2))
%    L(a)      = f(X,a) = product_by_n( f(Xn,a) )
%              = (2*pi*s^2)^(-N/2) * PI{1/Xn} * exp(-sum((log(Xn)-m)^2)/(2*s^2))
%    log(L(a)) = -N/2*log(2*pi) - N*log(s) - sum(log(Xn)) - sum((log(Xn)-m)^2)/(2*s^2)
%
%    The maximum likelihood point is found by the derivative of log(L(a)) with respect to "a":
%
%    diff(log(L(a)),s^2) = N/(2*s^4) * ( sum( (log(Xn)-m)^2 )/N - s^2 )
%                        = J(s^2) * (s^2_estimation - s^2)  
%    diff(log(L(a)),m)   = 2*N/(s^2) * ( sum(log(Xn))/N - m )
%                        = J(m) * (m_estimator - m)
%
%    Therefore, the (efficient) estimators are given by:
%
%               m   = sum( log( Xn ) )/N
%               s^2 = sum(( log(Xn)-m)^2)/N
%
%    The Cram?r-Rao Bounds for these estimator are:
%
%               VAR( m )   = 1/J(m)   = (s^2) / N
%               VAR( s^2 ) = 1/J(s^2) = (2*s^4) / N
%
%    NOTE: the ML estimator does not detect a deviation from the model.
%          therefore, check the RMS value !
%

if (nargin<1)
    error( 'fit_ML_log_normal - insufficient input arguments' );
end

% Estimation
% =============
lnx     = log(x(:));                 % should be column vectors !
N       = length(x);
m       = sum( lnx )/N;
s2      = (lnx-m)'*(lnx-m)/N;
CRB_m   = s2 / N;
CRB_s2  = (2*s2^2) / N;
[n,x_c] = hist( x,100 );
n       = n / sum(n*abs(x_c(2)-x_c(1)));
y       = sqrt(1/(2*pi))./(sqrt(s2)*x_c).*exp(- ((log(x_c)-m).^2)/(2*s2));
RMS     = sqrt( (y-n)*((y-n)')/ (x_c(2)-x_c(1))^2 / (length(x_c)-1) );

% finish summarizing results
% ============================
result = struct( 'm',m,'s2',s2,'CRB_m',CRB_m,'CRB_s2',CRB_s2,'RMS',RMS,'type','ML' );

% plot distribution if asked for
% ===============================
if (nargin>1)
    xspan = linspace(min(x),max(x),100);
    if ishandle( hAx )
        plot_log_normal( xspan,result,hAx,1 );
    else
        figure;
        plot_log_normal( xspan,result,gca,1 );
    end
end
