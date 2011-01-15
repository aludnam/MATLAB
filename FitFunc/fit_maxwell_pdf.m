function result = fit_maxwell_pdf( x,y,W,hAx )
% fit_maxwell_pdf - Non Linear Least Squares fit of the maxwellian distribution.
%                   given the samples of the histogram of the samples, finds the 
%                   distribution parameter that fits the histogram samples.
%
%    fits data to the probability of the form: 
%        p(r) = sqrt(2/pi)*(a^(-3/2))*(r^2)*exp(-(r^2)/(2*a))
%    with parameter: a
%
% format:   result = fit_maxwell_pdf( x,y,W,hAx )
%
% input:    y   - vector, samples of the histogram to be fitted
%           x   - vector, position of the samples of the histogram (i.e. y = f(x,a))
%           W   - matrix or scalar, a square weighting matrix of the size NxN where
%                 N = length(y), or 0 to indicate no weighting is needed.
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      a   - fitted parameter
%                      VAR - variance of the estimation
%                      type- weighted LS or not weighted LS
%                      iter- number of iteration for the solution
%

%
% Algorithm
% ===========
%
% We use the WLS algorithm to estimate the PDF from the samples.%
% The maxwell distribution is given by:
%
%    p(x,a) = sqrt(2/pi)*(a^(-3/2))*(x.^2).*exp(-(x.^2)/(2*a))
%           = Const * (a^(-3/2)) .* exp(-(x.^2)/(2*a))
%
%    note that X is known and therefore, considered a constant vector
%
% The non liner WLS estimator is given by:
%
%    a(n+1) = a(n) + inv(H'*W*H)*(H') * (y-h) = a(n) + G * err
% 
%    where:   h = p(x,a)
%             H = diff( p(x,a) ) with respect to "a"
%             W = weighting matrix of size NxN  (N = length(y))
%             a = a single parameter to be estimated
%
% The error estimation is given by:
%
%    VAR( a ) = G * VAR( err ) * (G')
%
%       or when W=I and the noise is a gaussian noise 
%
%    VAR( a ) = inv( H' * H )
%


if (nargin<3)
    error( 'fit_maxwell_pdf - insufficient input arguments' );
end

a       = x(find(y==max(y)))^2;         % initial guess
y       = y(:);                         % both should be column vectors !
x       = x(:);
x2      = x.^2;                         % save computation time
C       = sqrt(2/pi)*x2;                % a constant vector
thresh  = 0.995;                        % convergence threshold for the loop
last_cnt= inf;
iter    = 0;

% check weight matrix input
if (size(W,1)==length(y)) & (size(W,2)==length(y))
    weights_flag    = 1;
    type            = 'WLS';
else
    weights_flag    = 0;
    type            = 'LS';
end


% Estimation
% =============
if (weights_flag)
    % loop for convergence (with weighting matrix)
    % =============================================
    while (1)
        iter    = iter + 1;
        h       = C*(a^(-1.5)).*exp(-x2/(2*a));
        H       = h.*( x2/(2*a^2) - 3/(2*a) );
        HTW     = H'*W;
        e       = inv( HTW * H ) * HTW * (y-h);
        a       = a + e;
        control = e*e;
        if ( control > (last_cnt * thresh) )
            break;
        else
            last_cnt = control;
        end
    end

    % summarize results
    h           = C*(a^(-1.5)).*exp(-x2/(2*a));
    H           = h.*( x2/(2*a^2) - 3/(2*a) );
    HTW         = H'*W;
    G           = inv( HTW * H ) * HTW;
    err         = ( y - h );
    result.a    = a;
    result.VAR  = G * var( err ) * (G');
    result.RMS  = sqrt( (err')*err/ (x(2)-x(1))^2 / (length(err)-1) );
    result.iter = iter;
    result.type = type;
else

    % loop for convergence (without a weighting matrix) - assume white noise
    % ========================================================================
    while (1)
        iter    = iter + 1;
        h       = C*(a^(-1.5)).*exp(-x2/(2*a));
        H       = h.*( x2/(2*a^2) - 3/(2*a) );
        HT      = H';
        control = inv( HT * H );
        a       = a + control * HT * (y-h);
        if ( control>(last_cnt * thresh) )
            break;
        else
            last_cnt = control;
        end
    end
    
	% summarize results
    h           = C*(a^(-1.5)).*exp(-x2/(2*a));
	H           = h.*( x2/(2*a^2) - 3/(2*a) );
    err         = ( y - h );
    result.a    = a;
    result.VAR  = inv( (H') * H );
    result.RMS  = sqrt( (err')*err/ (x(2)-x(1))^2 / (length(err)-1) );
    result.iter = iter;
    result.type = type;
end


% plot distribution if asked for
% ===============================
if (nargin>3)
    if ishandle( hAx )
        plot_maxwell( x,result,hAx,2 );
    else
        figure;
        plot_maxwell( x,result,gca,2 );
    end
end
