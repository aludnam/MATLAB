function [out, c1, c2] = roiCornerCoord(in, topLeft, bottomRight, showim)
% Selects region of interest of the input image from the top left and bottom right coordinates.
% [out, c1, c2] = roiCornerCoord(in, topLeft, bottomRight, showim)
%
% in: image in
% topLeft: top-left coordinates of the ROI
% bottomRight: bottom-right coordinates of the ROI
% showim: show the original image with roi (default showim = 1)
% out: cut region of interest
% c1, c2: coordinates of the ROI corners

if ~exist('showim','var'); showim = 1; end
isdip = strcmp(class(in), 'dip_image'); 
if ~isdip
    in = dip_image(in); 
end
in = dip_image(in); 
dipshow(in); 
c1=topLeft;
c2=bottomRight;
out = in(c1(1):c2(1), c1(2):c2(2), :);
if ~isdip
    out = double(out);
end
if showim
    dipshow(in)
    line([c1(1) c1(1)], [c1(2) c2(2)])
    line([c2(1) c2(1)], [c1(2) c2(2)])
    line([c1(1), c2(1)], [c1(2) c1(2)])
    line([c1(1), c2(1)], [c2(2) c2(2)])
end
