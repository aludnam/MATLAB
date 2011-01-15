function h=vline(varargin)
%VLINE draws vertical lines
%   VLINE(X) draws one vertical line for each element in vector X
%   VLINE(AX,X) draws the lines to the axes specified in AX
%   VLINE(X,...) accepts HG param/value pairs for line objects
%   H=VLINE(...) returns a handle to each line created

%   Copyright 2007-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2008/05/23 15:39:00 $


if isempty(varargin) || (ishandle(varargin{1}) && length(varargin)==1)
    error('stats:hline:NotEnoughArgs','Not enough arguments');
end

if ishandle(varargin{1})
    ax=varargin{1};
    varargin=varargin(2:end);
else
    ax=gca;
end

x = varargin{1};
varargin=varargin(2:end);
    
hh=graph2d.constantline(x,'parent',ax,varargin{:});
hh = changedependvar(hh,'x');

if nargout>0
    h=hh;
end