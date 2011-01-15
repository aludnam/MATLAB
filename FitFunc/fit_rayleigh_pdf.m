function result = fit_rayleigh_pdf( x,y,W,hAx )
% fit_rayleigh_pdf - Non Linear Least Squares fit of the Rayleigh distribution.
%                    given the samples of the histogram of the samples, finds the 
%                    distribution parameter that fits the histogram samples.
%
%    fits data to the probability of the form: 
%        p(r)=r*exp(-r^2/(2*s))/s
%    with parameter: s
%
% format:   result = fit_rayleigh_pdf( x,y,W,hAx )
%
% input:    y   - vector, samples of the histogram to be fitted
%           x   - vector, position of the samples of the histogram (i.e. y = f(x,a))
%           W   - matrix or scalar, a square weighting matrix of the size NxN where
%                 N = length(y), or 0 to indicate no weighting is needed.
%           hAx - handle of an axis, on which the fitted distribution is plotted
%                 if h is given empty, a figure is created.
%
% output:   result  - structure with the fields
%                      s   - fitted parameter
%                      VAR - variance of the estimation
%                      type- weighted LS or not weighted LS
%                      iter- number of iteration for the solution
%

%
% Algorithm
% ===========
%
% We use the WLS algorithm to estimate the PDF from the samples.%
% The rayleigh distribution is given by:
%
%    p(x,s) = x*exp(-x^2/(2*s))/s
%
%    note that X is known and therefore, considered a constant vector
%
% The non liner WLS estimator is given by:
%
%    s(n+1) = s(n) + inv(H'*W*H)*(H') * (y-h) = s(n) + G * err
% 
%    where:   h = p(x,s)
%             H = diff( p(x,a) ) with respect to "s"
%             W = weighting matrix of size NxN  (N = length(y))
%             s = a single parameter to be estimated
%
% The error estimation is given by:
%
%    VAR( s ) = G * VAR( err ) * (G')
%
%       or when W=I and the noise is a gaussian noise 
%
%    VAR( s ) = inv( H' * H )
%

if (nargin<3)
    error( 'fit_rayleigh_pdf - insufficient input arguments' );
end

s       = x(find(y==max(y))).^2;        % initial guess
y       = y(:);                         % both should be column vectors !
x       = x(:);
x2      = x.^2;                         % save computation time
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
        h       = x.*exp( -x2/(2*s) )/s;
        H       = h.*(x2/(2*s^2) - 1/s);
        HTW     = H'*W;
        e       = inv( HTW * H ) * HTW * (y-h);
        s       = s + e;
        control = e*e;
        if ( control > (last_cnt * thresh) )
            break;
        else
            last_cnt = control;
        end
    end

    % summarize results
    h       = x.*exp( -x2/(2*s) )/s;
    H       = h.*(x2/(2*s^2) - 1/s);
    HTW     = H'*W;
    G       = inv( HTW * H ) * HTW;
    err     = (y-h);
    VAR     = G * var(err) * (G');
else

    % loop for convergence (without a weighting matrix)
    % ==================================================
    while (1)
        iter    = iter + 1;
        h       = x.*exp( -x2/(2*s) )/s;
        H       = h.*(x2/(2*s^2) - 1/s);
        HT      = H';
        control = inv( HT * H );
        s       = s + control * HT * (y-h);
        if ( control>(last_cnt * thresh) )
            break;
        else
            last_cnt = control;
        end
    end
    
	% summarize results
    h       = x.*exp( -x2/(2*s) )/s;
    H       = h.*(x2/(2*s^2) - 1/s);
    err     = (y-h);
    VAR     = inv( (H') * H );
end


% finish summarizing results
% ============================
result.s    = s;
result.VAR  = VAR;
result.RMS  = sqrt( (err')*err/ (x(2)-x(1))^2 / (length(err)-1) );
result.iter = iter;
result.type = type;

% plot distribution if asked for
% ===============================
if (nargin>3)
    if ishandle( hAx )
        plot_rayleigh( x,result,hAx,4 );
    else
        figure;
        plot_rayleigh( x,result1,gca,4 );
    end
end

