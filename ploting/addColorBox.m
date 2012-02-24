function addColorBox(topLeft,bottomRight,color,linewidth,flipCoord)
% Adds color box to an image. 
% addColorBox(topLeft,bottomRight,color,linewidth,flipCoord)
% topLeft, bottomRight are the coorinated of the box corners. 
% color: color of hte box. default: color = 'red'
% linewidth: width of the line. default: linewidth=1
% flipCoord: if st to 1 then coordinates are fliped (this is for marix in double)

if ~exist('color','var'); color = 'red'; end
if ~exist('linewidth','var'); linewidth= 1; end
if ~exist('flipCoord','var');flipCoord=0; end


c1=topLeft;
c2=bottomRight;
if flipCoord % this is because of hte dip_image mismatch of coordiantes system
    c1=fliplr(c1);
    c2=fliplr(c2); 
end

line([c1(1) c1(1)], [c1(2) c2(2)],'color',color,'linewidth',linewidth)
line([c2(1) c2(1)], [c1(2) c2(2)],'color',color,'linewidth',linewidth)
line([c1(1), c2(1)], [c1(2) c1(2)],'color',color,'linewidth',linewidth)
line([c1(1), c2(1)], [c2(2) c2(2)],'color',color,'linewidth',linewidth)
