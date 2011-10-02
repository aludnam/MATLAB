function [Hscale,Htext]=plot_scale(Pos,Scale,Length,Color,UnitsName,Orient);
%-------------------------------------------------------------------
% I slightly modified this file so that the widht of the bar
% (...,'linewidth',4,...) and prints no text at the bar
% plot_scale function     Add a scale bar on a plot or image.
% Input  : - Scale bar position [X, Y].
%          - Scale [Units per pixel].
%          - Scale bar length in units (e.g., arcsec).
%          - Color, default is 'k';
%          - Units name, default is 'arcsec'.
%          - Scale bar orientation:
%            'h' - Horizontal (default).
%            'v' - Vertical. 
% Output : - Handle for the scale line.
%          - Handle for the text.
% Tested : Matlab 7.0
%     By : Eran O. Ofek          July 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-------------------------------------------------------------------
DistFactor = 0.05;
if (nargin==3),
   Color     = 'k';
   UnitsName = 'arcsec';
   Orient    = 'h';
elseif (nargin==4),
   UnitsName = 'arcsec';
   Orient    = 'h';
elseif (nargin==5),
   Orient    = 'h';
elseif (nargin==6),
   % do nothing
else
   error('Illegal Number of input arguments');
end

NextPlot = get(gca,'NextPlot');
hold on;

XLim   = get(gca,'XLim');
YLim   = get(gca,'YLim');
Xdiff  = abs(diff(XLim));
Ydiff  = abs(diff(YLim));

switch Orient
 case 'h'
    LineX = Pos(1) + 0.5.*Length./Scale.*[-1;+1];
    LineY = Pos(2).*[+1;+1];
    DistXdir = 0;
    DistYdir = -1;
 case 'v'
    LineX = Pos(1).*[+1;+1];
    LineY = Pos(2) + 0.5.*Length./Scale.*[-1;+1];
    DistXdir = +1;
    DistYdir = 0;
 otherwise
    error('Unknown Orient Option');
end

%--- plot line ---
Hscale = plot(LineX,LineY);
set(Hscale,'Color',Color,'linewidth',4);

%--- plot text ---
DistX    = DistXdir.*DistFactor.*Xdiff;
DistY    = DistYdir.*DistFactor.*Ydiff;

% Htext    = text(Pos(1)+DistX,Pos(2)+DistY,sprintf('%5.1f %s',Length,UnitsName));
Htext    = text(Pos(1)+DistX,Pos(2)+DistY,sprintf('%5.1f %s',[],[]));

switch Orient,
 case 'h'
    % do nothing
 case 'v'
    set(Htext,'Rotation',90);
 otherwise
    error('Unknown Orient Option');
end
set(Htext,'HorizontalAlignment','center','Color',Color);

set(gca,'XLim',XLim);
set(gca,'YLim',YLim);

set(gca,'NextPlot',NextPlot);