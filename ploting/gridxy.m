function hh = gridxy(x,varargin)
% GRIDXY - Plot grid lines
%   GRIDXY(X) plots vertical grid lines at the positions specified
%   by X. GRIDXY(X,Y) also plots horizontal grid lines at the positions
%   specified by Y. GRIDXY uses the current axes, if any. Lines outside
%   the plot area are plotted but not shown. When X or Y is empty no vertical
%   or horizontal lines are plotted.
%
%   The lines are plotted as a single graphics object. H = GRIDXY(..) returns
%   a graphics handle to that line object. 
%
%   GRIDXY(..., 'Prop1','Val1','Prop2','Val2', ...) uses the properties
%   and values specified for color, linestyle, etc. Execute GET(H), where H is
%   a line handle, to see a list of line object properties and their current values.
%   Execute SET(H) to see a list of line object properties and legal property values.
%
%   Examples
%     % some random plot
%       plot(10*rand(100,1), 10*rand(100,1),'bo') ; 
%     % horizontal red dashed grid
%       gridxy([1.1 3.2 4.5],'Color','r','Linestyle',':') ;
%     % vertical solid thicker yellowish grid, and store the handle
%       h = gridxy([],[2.1:0.7:5 8],'Color',[0.9 1.0 0.2],'linewidth',3) ;
%
%   GRIDXY can be used to plot a irregular grid on the axes.
%
%   See also PLOT, REFLINE, GRID, AXES, REFLINEXY

% NOTE: This function was previously known as XYREFLINE

% for Matlab R13
% version 2.2 (feb 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% Created (1.0) feb 2006
% 2.0 apr 2007 - renamed from reflinexy to gridxy, reflinexy is now used
%               for plotting intersection between X and Y axes
% 2.1 apr 2007 - add error check for line properties
% 2.2 feb 2008 - added set(gca,'layer','top') to put gridlines behind the
%                axis tick marks

error(nargchk(1,Inf,nargin)) ;

% check the arguments
if ~isnumeric(x),
    error('Numeric argument expected') ;
end

if nargin==1,
    y = [] ;
    va = [] ;
else
    va = varargin ;
    if ischar(va{1}),
        % optional arguments are
        y = [] ;
    elseif isnumeric(va{1})        
        y = va{1} ;
        va = va(2:end) ;
    else
        error('Invalid second argument') ;
    end
    if mod(size(va),2) == 1,
        error('Property-Value have to be pairs') ;
    end
end

% get the axes to plot in
hca=get(get(0,'currentfigure'),'currentaxes');
if isempty(hca),
    warning('No current axes found') ;
    return ;
end

% get the current limits of the axis
% used for limit restoration later on
xlim = get(hca,'xlim') ;
ylim = get(hca,'ylim') ;

% setup data for the vertical lines
xx1 = repmat(x(:).',3,1) ;
yy1 = repmat([ylim(:) ; nan],1,numel(x)) ;

% setup data for the horizontal lines
xx2 = repmat([xlim(:) ; nan],1,numel(y)) ;
yy2 = repmat(y(:).',3,1) ;


% create data for a single line object
xx1 = [xx1 xx2] ;
if ~isempty(xx1),     
    yy1 = [yy1 yy2] ;
    % add the line to the current axes
    np = get(hca,'nextplot') ;
    set(hca,'nextplot','add') ;
    h = line('xdata',xx1(:),'ydata',yy1(:)) ;          
    set(hca,'ylim',ylim,'xlim',xlim) ; % reset the limits
    
    uistack(h,'bottom') ; % push lines to the bottom of the graph
    set(hca,'nextplot',np,'Layer','top') ;    % reset the nextplot state

    if ~isempty(va),
        try
            set(h,va{:}) ; % set line properties        
        catch
            msg = lasterror ;
            error(msg.message(21:end)) ;
        end
    end

else
    h = [] ;
end

if nargout==1,     % if requested return handle
    hh = h ;
end









