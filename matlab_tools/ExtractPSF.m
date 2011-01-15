im = readtimeseries('beads1_000');
prj=squeeze(sum(im(:,:,3:5),[],3))

psf=extract(prj,[17 17],[181,91])
border=(abs(xx(size(psf))) >= size(psf,1)/2-1) | (abs(yy(size(psf))) >= size(psf,2)/2-1)
psf=psf-mean(psf,border);
psf=DampEdge(psf,0.12,2,1,2)

writeim(psf,'psf1');

