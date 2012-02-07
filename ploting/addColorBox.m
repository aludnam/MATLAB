function addColorBox(topLeft,bottomRight,color,linewidth)
% Adds color box to an image. 
% addColorBox(topLeft,bottomRight,color,linewidth)
% topLeft, bottomRight are the coorinated of the box corners. 
% color: color of hte box. default: color = 'red'
% linewidth: width of the line. default: linewidth=1

if ~exist('color','var'); color = 'red'; end
if ~exist('linewidth','var'); linewidth= 1; end
c1=topLeft;
c2=bottomRight;

line([c1(1) c1(1)], [c1(2) c2(2)],'color',color,'linewidth',linewidth)
line([c2(1) c2(1)], [c1(2) c2(2)],'color',color,'linewidth',linewidth)
line([c1(1), c2(1)], [c1(2) c1(2)],'color',color,'linewidth',linewidth)
line([c1(1), c2(1)], [c2(2) c2(2)],'color',color,'linewidth',linewidth)
