function [x,y]=sample2dsubpix(image, n)
% [x,y]=sample2dsubpix(image, n)
% samples from the distribution provided as a 2d pixelised intensity image
% image - pixelised intensity image
% n - number of samples
if ~exist('n', 'var')
    n=1;
end
[xpix,ypix]=sample2d(image,n);
x=xpix-1 + rand;
y=ypix-1 + rand;