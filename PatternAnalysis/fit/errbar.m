function errorbar(x, y, l,u)
%ERRBAR Plot error bars.
%	ERRBAR(X,Y,L,U) plots errorbars on an existing graph of vector X 
%	vs. vector Y with error bars specified by the vectors L and U.  The vectors 
%	X,Y, L and U must be the same length.  If X,Y, L and U are matrices then 
%	each column produces a separate line.  The error bars are each drawn a 
%	distance of U(i) above and L(i) below the points in (X,Y) so that each bar 
%	is L(i) + U(i) long.  
%	ERRBAR(X,Y,L) puts error bars on a graph of X vs Y that are symmetric  
%	about Y and are 2*L(i) long.
%
%	ERRBAR(Y,L) puts error bars [Y-L Y+L] on a graph of Y.
%
%	For example,
%
%	   x = 1:10;
%	   y = sin(x);
%	   e = std(y)*ones(size(x));
%	   plot(x,y);
%	   errbar(x,y,e)
%
%	draws symmetric error bars of unit standard deviation.

%	L. Shure 5-17-88, 10-1-91 B.A. Jones 4-5-93
%	Copyright (c) 1984-94 by The MathWorks, Inc.

hold on
if min(size(x))==1,
  npt = length(x);
  x = x(:);
  y = y(:);
  if nargin == 3,  
			l = l(:);
		elseif nargin == 4 | nargin == 5
			l = l(:);
			u = u(:);
		end
else
  [npt,n] = size(x);
end

if nargin == 3,  
	u = l;
end

if nargin == 2
  l = y;
		u = y;
  y = x;
  [m,n] = size(y);
  x(:) = [1:npt]'*ones(1,n);;
end
if isstr(x) | isstr(y) | isstr(l) | isstr(u)
	error('Arguments must be numeric.')
end

if any(size(x)~=size(y)) | any(size(x)~=size(l)) |  any(size(x)~=size(u)),
  error('The sizes of X, Y, L and U must be the same.');
end

tee = (max(x(:))-min(x(:)))/100;  % make tee .02 x-distance for error bars
xl = x - tee;
xr = x + tee;
n = size(y,2);

% Plot graph and bars
cax = newplot;
next = lower(get(cax,'NextPlot'));
% build up nan-separated vector for bars
xb = [];
yb = [];
nnan = nan*ones(1,n);
for i = 1:npt
    ytop = y(i,:) + u(i,:);
    ybot = y(i,:) - l(i,:);
    xb = [xb; x(i,:); x(i,:) ; nnan; xl(i,:);xr(i,:);nnan ;xl(i,:);xr(i,:);nnan];
    yb = [yb; ytop;ybot;nnan;ytop;ytop;nnan;ybot;ybot;nnan];
end

plot(xb,yb)
hold off
