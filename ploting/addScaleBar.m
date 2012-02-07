function addScaleBar(in,barLength,barWidth, barPosition,barColor)
% Adds scale bar into an image. 
% out=addScaleBar(in,barLength,barPosition,barColor)
% in: image to which the scale bar should be added. 
% barLength: length of the scale bar in pixels (defaults barLength=10);
% barWidth: thickness of the scale bar in pixesl (default: barWidth = 3); 
% barPosition: position of the scale bar in the image (default: barPosition = size(in)*9/10;
% barColor: color of the scale bar (default: 'white')

if and(exist('in','var'),~isempty(in)); 
    dipshow(in); 
else
    % added to a current picture
    c=get(gca);
    in = zeros(floor([c.XLim(2),c.YLim(2)])); % just to get size right for barPosition
end;
if ~exist('barLength','var'); barLength = 10; end; 
if ~exist('barWidth','var'); barWidth = 3; end;
if ~exist('barPosition','var'); barPosition = size(in)*9/10; end;
if ~exist('barColor','var'); barColor = [1 1 .99]; end;



line([barPosition(1) barPosition(1)-barLength], [barPosition(2) barPosition(2)],'color',barColor,'linewidth',barWidth)

    