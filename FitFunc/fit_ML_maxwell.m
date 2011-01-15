function result = fit_ML_maxwell( x,hAx )
% fit_ML_maxwell - Maximum Likelihood fit of the maxwellian distribution of i.i.d. samples!.
%                  Given the samples of a maxwellian distribution, the PDF parameter is found
%
%    fits data to the probability of the form: 
%        p(r) = sqrt(2/pi)*(a^(-3/2))*(r^2)*exp(-(r^2)/(2*a))
%    with parameter: a
%
% format:   result = fit_ML_maxwell( x,hAx )
%
% input:    x   - vector, samples with maxwellian distribution to be parameterized
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      a   - fitted parameter
%                      CRB - Cram?r-Rao Bound for the estimator value
%                      RMS - RMS error of the estimation 
%                      type- 'ML'
%

%
% Algorithm
% ===========
%
% We use the ML algorithm to estimate the PDF from the samples.
% The maxwell destribution is given by:
%
%    p(x,a) = sqrt(2/pi)*(a^(-3/2))*(x.^2).*exp(-(x.^2)/(2*a))
%
%    where x are the samples which distribute by the function p(x,a)
%            and are assumed to be i.i.d !!!
%
% The ML estimator is given by:
%
%    f(Xn,a)   = sqrt(2/pi)*a^(-3/2)*(Xn^2)*exp( -(Xn^2)/(2*a) )
%    L(a)      = f(X,a) = product_by_n( f(Xn,a) )
%              = (2/pi)^(N/2) * a^(-3*N/2) * PI((Xn^2)) * exp( -sum(Xn^2)/(2*a) )
%    log(L(a)) = N/2*log(2/pi) - 3*N/2*log(a) + 2*sum(log(Xn)) - sum(Xn^2)/(2*a)
%
%    The maximum likelihood point is found by the derivative of log(L(a)) with respect to "a":
%
%    diff(log(L(a))) = sum(Xn^2)/(2*a^2) - 3*N/(2*a) = (3*N)/(2*a^2) * ( sum(Xn^2)/(3*N) - a )
%                    = J(a) * (a_estimation - a)  
%
%    Therefore, the (efficient) estimator is given by:
%
%               a = sum( Xn^2 ) / (3 * N)
%
%    The Cram?r-Rao Bound for this estimation is:
%
%               VAR( a ) = 1/J(a) = (2*a^2)/(3*N)
%
%    NOTE: the ML estimator does not detect a deviation from the model.
%          therefore, check the RMS value !
%

if (nargin<1)
    error( 'fit_ML_maxwell - insufficient input arguments' );
end

% Estimation
% =============
x       = x(:);                 % should be column vectors !
N       = length(x);
a       = sum(x.^2)/(3*N);
CRB     = (2*a^2)/(3*N);
[n,x_c] = hist( x,100 );
n       = n / sum(n*abs(x_c(2)-x_c(1)));
y       = sqrt(2/pi)*(a^(-3/2))*(x_c.^2).*exp(-(x_c.^2)/(2*a));
RMS     = sqrt( (y-n)*((y-n)')/ (x_c(2)-x_c(1))^2 / (length(x_c)-1) );

% finish summarizing results
% ============================
result = struct( 'a',a,'CRB',CRB,'RMS',RMS,'type','ML' );

% plot distribution if asked for
% ===============================
if (nargin>1)
    xspan = linspace(min(x),max(x),100);
    if ishandle( hAx )
        plot_maxwell( xspan,result,hAx,1 );
    else
        figure;
        plot_maxwell( xspan,result,gca,1 );
    end
end
