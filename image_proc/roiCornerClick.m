function [out, c1, c2] = roiCornerClick(in, showim, color)
% Selects region of interest of the input image by clicking on the ROI
% corners.
% [out, c1, c2] = roiCornerClick(in, showim, color)
%
% in: image in
% showim: show the original image with roi (default showim = 1)
% color: [string] color of the line (default color = 'blue') 
% out: cut region of interest
% c1, c2: coordinates of the ROI corners

if ~exist('showim','var'); showim = 1; end
if ~exist('color','var'); color='blue'; end
isdip = strcmp(class(in), 'dip_image'); 
if ~isdip
    in = dip_image(in); 
end
dipshow(in); 
fprintf('Click to upper left corner.\n');
c1=dipgetcoords;
fprintf('Click to bottom right corner.\n');
c2=dipgetcoords;
c1=c1(1:2);
c2=c2(1:2);
if ndims(in)>2
    out = in(c1(1):c2(1), c1(2):c2(2), :);
else 
    out = in(c1(1):c2(1), c1(2):c2(2));
end
fprintf('Selected UPPER LEFT    coordinates are:[%g, %g]\n',c1)
fprintf('Selected BOTTM RIHGT   coordinates are:[%g %g]\n',c2)
    
if ~isdip
    out = double(out);
end
if showim
    dipshow(in)
    addColorBox(c1,c2,color);   
end
