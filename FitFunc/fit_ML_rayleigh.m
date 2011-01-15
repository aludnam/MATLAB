function result = fit_ML_rayleigh( x,hAx )
% fit_ML_rayleigh - Maximum Likelihood fit of the rayleigh distribution of i.i.d. samples!.
%                  Given the samples of a rayleigh distribution, the PDF parameter is found
%
%    fits data to the probability of the form: 
%        p(r)=r*exp(-r^2/(2*s))/s
%    with parameter: s
%
% format:   result = fit_ML_rayleigh( x,hAx )
%
% input:    x   - vector, samples with rayleigh distribution to be parameterized
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      s   - fitted parameter
%                      CRB - Cram?r-Rao Bound for the estimator value
%                      RMS - RMS error of the estimation 
%                      type- 'ML'
%

%
% Algorithm
% ===========
%
% We use the ML algorithm to estimate the PDF from the samples.
% The rayleigh destribution is given by:
%
%    p(x,s) = x / s * exp(-x^2/(2*s))
%
%    where x are the samples which distribute by the function p(x,s)
%            and are assumed to be i.i.d !!!
%
% The ML estimator is given by:
%
%    f(Xn,s)   = Xn / s * exp(-Xn^2/(2*s))
%    L(s)      = f(X,s) = product_by_n( f(Xn,s) )
%              = PI(Xn) * (s^(-N)) * exp( -sum(Xn^2)/(2*s) )
%    log(L(s)) = sum(log(Xn)) - N*log(s) - sum(Xn^2)/(2*s)
%     
%    The maximum likelihood point is found by the derivative of log(L(s)) with respect to "s":
%
%    diff(log(L(s))) = -N/s + sum(Xn^2)/(2*s^2) = (N/s^2) * ( sum(Xn^2)/(2*N) - s ) = 
%                    = J(s) * (s_estimation - s)  
%
%    Therefore, the (efficient) estimator is given by:
%
%               s = sum( Xn^2 ) / (2 * N)
%
%    The Cram?r-Rao Bound for this estimation is:
%
%               VAR( s ) = 1/J(s) = (s^2)/N
%
%    NOTE: the ML estimator does not detect a deviation from the model.
%          therefore, check the RMS value !
%

if (nargin<1)
    error( 'fit_ML_rayleigh - insufficient input arguments' );
end

% Estimation
% =============
x       = real(x(:));                 % should be column vectors !
N       = length(x);
s       = sum(x.^2)/(2*N);
CRB     = (s^2)/N;
[n,x_c] = hist( x,100 );
n       = n / sum(n*abs(x_c(2)-x_c(1)));
y       = x_c.*exp(-x_c.^2/(2*s))/s;
RMS     = sqrt( (y-n)*((y-n)')/ (x_c(2)-x_c(1))^2 / (length(x_c)-1) );

% finish summarizing results
% ============================
result = struct( 's',s,'CRB',CRB,'RMS',RMS,'type','ML' );

% plot distribution if asked for
% ===============================
if (nargin>1)
    xspan = linspace(min(x),max(x),100);
    if ishandle( hAx )
        plot_rayleigh( xspan,result,hAx,3 );
    else
        figure;
        plot_rayleigh( xspan,result,gca,3 );
    end
end
