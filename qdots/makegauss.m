function wg = makegauss(c, sig, sizevec)
% function wg = makegauss(c, sig, sizevec)

nx=sizevec(1);
ny=sizevec(2);

ng = length(c)/2;
[X,Y] = meshgrid(0:nx-1,0:ny-1);
n0=1/(2*pi*sig^2);
cr = reshape(c,2,ng)';

for ii=1:ng
    wg(:,:,ii)=n0*exp(-((X-cr(ii,1)).^2+(Y-cr(ii,2)).^2)/(2*sig^2));
end