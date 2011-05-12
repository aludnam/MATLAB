function [imout, xrout, yrout] = convolvepoints(x,y,psf_tmp,np,pixelsize)
% [imout, xrout, yrout] = convolvepoints(x,y,psf_tmp,np,pixelsize)
%
% Generates an image with multiple (N) PSFs placed at coordinates ([x,y]).
% x,y:          x,y coordiantes (each Nx1 vector)
% psf_tmp:      image of a psf
% np:           vector of intensities of each PSF (default: 1 for all)
% pixelsize:    size of a pixel if x and y are in [pixels], if x and
% y are in nm then pixelsize = 1; (default: pixelsize = 1)
%
% Size of the psf image should be odd with the maximum in the central pixel. 

if ~exist ('np', 'var')
    np = ones(size(x));
end

if ~isequal(length(x), length(np))
    np = np * ones(size(x));
end

if ~exist ('pixelsize', 'var')
    pixelsize = 1;
end

sx = size(x);
sy = size(y);
if sx(1)<sx(2); x =x'; end %to make sure that x and y are column vectors
if sy(1)<sy(2); y =y'; end

x = x/pixelsize; % all in pixels...
y = y/pixelsize; 
minx = min(x); maxx = max(x);
miny = min(y); maxy = max(y);
margin = [5 5]; % margin around the image - provides spave for shifting to the very edge

psf_tmpnorm = double(psf_tmp/sum(psf_tmp(:))); %normalize to 1
sizeim_tmp = ceil([maxx-minx, maxy-miny])+2*margin;
margincorr = 1-mod(sizeim_tmp,2); % correction of the margin to make it odd...
sizeim = sizeim_tmp + margincorr; 

centerpoints = [maxx-minx, maxy-miny]/2;
centerim = sizeim/2;
psf = padimage(psf_tmpnorm,sizeim);
xr=x-minx+.5; % so that it starts from .5 (point has center in the middle of hte first pixel)
yr=y-miny+.5;

imout_tmp = zeros(size(psf));
for ii=1:length(x)
    imout_tmp = imout_tmp + np(ii)*FourierShift2D(psf, [xr(ii), yr(ii)]-centerpoints);
end
imout = imout_tmp(margin(1):end-(margin(1)+margincorr(1)),margin(2):end-(margin(2)+margincorr(2)));
pointscorr=centerim-centerpoints-margin;
coorout = bsxfun(@plus, [xr, yr],pointscorr)+.5+1;
xrout = coorout(:,2); % some confusion in x and y (colum of the matrix are hte x coordiantes, rows are the y coordinates...)
yrout = coorout(:,1);