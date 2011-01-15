function plot_mix_gaussian( u,sig,prob,X )
%
% plot_mix_gaussian - plot the samples and the estimation. 
%                     for 1D distribution, plot the normalized histogram and the 1D distribution function        
%                     for 2D distribution, plot the samples and the contour of FWHM
%                     for mD distribution (m>2) - do not plot nothing.
%
% format:   plot_mix_gaussian( u,sig,prob,X )
%
% input:    u       - mean of each gaussian in the distribution. 1xM vector or 2xM matrix.
%                     each gaussian mean is stored in a separate column.
%           sig     - for 1D distribution -> standard deviation of each gaussian in the distribution
%                     each gaussian mean is stored in a separate column -> a 1xM vector
%                     for 2D distribution -> covariance matrix for each gaussian in the distribution
%                     2x2xM matrix, the 3rd dimension is the gaussians index,
%                     the 1st and 2nd dimensions are the covariance matrix
%           prob    - probability of each gaussian in the distribution to create the current sample.
%                     this is a 1xM vector
%           X       - the samples, 1xN vector or 2xN matrix, depends on the dimension of the 
%                     distribution (i.e. 1D or 2D)
%
% output:   to the graphic current axis.
%
%
%

% check input
if (nargin<3)
    error( 'plot_mix_gaussian - insufficient input parameters' );
end
if ~exist( 'X' )
    X = [];
end

% constants
nbins = 200;

% check the size of the input to determin if it's 1D or 2D distribution
if (size(u,1)==2) & (size(u,2)==size(prob,2))
    plot_mix_gaussian_2d( u,sig,prob,X );
else
    plot_mix_gaussian_1d( u,sig,prob,X,nbins );
end
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              Inner function implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_mix_gaussian_1d( u,sig,prob,X,nbins )
% plot normalized histogram and normalized distribution on top of it

% constant
points = 1000;

if ~isempty( X )
    [n,x]   = hist( X,nbins );              % calc the histogram
    dx      = x(2)-x(1);                    % calc a single bin width
    n       = n / sum( n*dx );              % normalize histogram to have area of 1
    bar( x,n,'hist' );                      % plot normalized histogram
    xlim( [x(1)-dx/2,x(end)+dx/2] );        % make sure that the axis is squeezed to it's limits
    g       = gca;
    c       = [1 0 0];                      % choose the red color
    x_c     = linspace( x(1)-dx/2,x(end)+dx/2,points );
else
    g       = gca;
    X       = xlim( g );                    % get the axis limits
    c       = [0 1 0];                      % choose the green color
    x_c     = linspace( X(1),X(2),points );
end

% plot the distribution
for m = 1:length(prob)
    y   = prob(m) / sqrt( 2*pi*sig(m)^2 ) * exp( -((x_c-u(m)).^2)/(2*sig(m)^2) );
    line( x_c,y,'color',c,'parent',g,'linewidth',2 );
end
shg;
drawnow;

% ----------------------------------------------------------------------------------------

function plot_mix_gaussian_2d( u,covar,prob,X )
% plot 2D samples and distribution data on top of it

if ~isempty( X )
    if size(X,2)<size(X,1),
        X = X.';
    end
    plot( X(1,:),X(2,:),'.k' );             % plot samples
    xlim( [min(X(1,:)) max(X(1,:))] );      % squeeze axis limits
    ylim( [min(X(2,:)) max(X(2,:))] );       % squeeze axis limits
    g       = gca;
    c       = [1 0 0];                      % choose the red color
else
    g       = gca;
    c       = [0 1 0];                      % choose the green color
end

% plot each gaussian's FWHM and mean (center)
a           = linspace(0,2*pi,361); % degree vector
unit_circle = [cos(a);sin(a)];      % a unit circle for the case of sig_x=sig_y=1
for m = 1:length(prob)
    U       = u(:,m);               % the mean of the mth gaussian
    [V,D]   = eig(covar(:,:,m));    % eigen value and vectors for the covariance of the mth gaussian
    A       = V*sqrt(D);            % factorization: covar = A*A' = positive matrix
    outline = A * unit_circle * sqrt(2*log(2)) + U*ones(1,length(unit_circle)); % the contour
    line( U(1),U(2),'color',c,'parent',g,'marker','o' );
    line( outline(1,:),outline(2,:),'color',c,'parent',g,'linewidth',2 );
end