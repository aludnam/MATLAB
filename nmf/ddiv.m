function f = ddiv(cr, varargin)
% gf = gradientd(c, varargin)
% c - position of each ceneter (2*Ncomp x 1 dimensional row vector)
% Vxt = varargin{1};  %data
% Hkt = varargin{2};  %estimated H from NMF update
% sig = varargin{3};  %variance of the psf gaussian approximation
% [nx, ny] = varargin{4}; %size of the image

Vxt = varargin{1};  %data
Hkt = varargin{2};  %estimated H from NMF update
sig = varargin{3};  %variance of the psf gaussian approximation
sizexy = varargin{4}; %size of the image
nx = sizexy(1);
ny = sizexy(2);

nc = length(cr)/2;  %number of components
sH = size(Hkt,1);

if nc<sH   %one component wac background....
    bg = Hkt(sH,1)/(nx*ny);
    Hkt = Hkt(1:sH-1,:); %only components, no background...
else
    bg = 0;
end

c=reshape(cr',2,nc)'; %nc x 2 dim -> centers of components

[X, Y] = meshgrid(0:nx-1, 0:ny-1);
Xnc = repmat(X,[1 1 nc]);
Ync = repmat(Y,[1 1 nc]);
cx = shiftdim(repmat(c(:,1), [1, ny, nx]),1);
cy = shiftdim(repmat(c(:,2), [1, ny, nx]),1);
Xc = Xnc-cx;
Yc = Ync-cy;

n0=1/(2*pi*sig^2);
Wxk_mat = n0*exp(-(Xc.^2+Yc.^2)/(2*sig^2));    %ny x nx x nc matrix - each slice is individual component...
Wxk = (reshape(Wxk_mat, nx*ny, nc));

fmat = Vxt.*log(Vxt./(Wxk*Hkt+bg)) - Vxt + (Wxk*Hkt+bg);
f = sum(fmat(:));

