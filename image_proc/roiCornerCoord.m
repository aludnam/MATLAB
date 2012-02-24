function [out, c1, c2] = roiCornerCoord(in, topLeft, bottomRight, showim, color, flipCoord)
% Selects region of interest of the input image from the top left and bottom right coordinates.
% [out, c1, c2] = roiCornerCoord(in, topLeft, bottomRight, showim, color, flipCoord)
%
% in: image in
% topLeft: top-left coordinates of the ROI
% bottomRight: bottom-right coordinates of the ROI
% showim: show the original image with roi (default showim = 1)
% color: [string] color of the line (default color = 'blue') 
% flipCoord: set to 1 for doubles
% out: cut region of interest
% c1, c2: coordinates of the ROI corners


if ~exist('showim','var'); showim = 1; end
if ~exist('color','var'); color='blue'; end
if ~exist('flipCoord','var');flipCoord=0; end

isdip = strcmp(class(in), 'dip_image'); 
if ~isdip
    in = dip_image(in); 
end
if flipCoord % this is because of hte dip_image mismatch of coordiantes system
    topLeft=fliplr(topLeft)-1;
    bottomRight=fliplr(bottomRight)-1; 
end

in = dip_image(in); 
c1=topLeft;
c2=bottomRight;
if ndims(in)>2
    out = in(c1(1):c2(1), c1(2):c2(2), :);
else 
    out = in(c1(1):c2(1), c1(2):c2(2));
end

if ~isdip
    out = double(out);
end
if showim
    dipshow(in)
    addColorBox(c1,c2,color)    
end
