function [mwcn,m]=analyzeWconfPSF(w,peval,p)
% [mwcn,m]=analyzeWconfPSF(w,peval,p)
% Finds local maxima in the results w (nx*ny,ncomp) convoluted with PSF
% (generated from parameters p.lambda, p.NA, p.ri,p.pixelsize)
% m = binary images of the local maxima. 
% mwcn = image of local maxima (non zero pixels) with values of hte convoluted image, nornalised to the maximum of each column.  

o = kSimPSF( {'lambdaEm',p.lambda;'Pi4Em',0;'relEmInt',1;'relEmPhase',0;'na',p.NA;'ri',p.ri;'sX',p.nx;'sY',p.ny;'sZ',1;'scaleX',p.pixelsize;'scaleY',p.pixelsize;'scaleZ',p.pixelsizeZ;'lambdaEx',488;'twophoton',0;'confocal',0;'nonorm',0;'pinhole',1;'Pi4Ex',0;'relExInt',1});
psf=double(o); % it is normalized to 1
wpix=reshape(w,peval.nx, peval.ny, peval.ncomp); % each frame normalized to 1
wpixc=(convstack(wpix,psf,'same'));
wc=reshape(wpixc,peval.nx*peval.ny,peval.ncomp);
mpix=maximastack(wpixc); % binary image of local maxima locations in each frame
m=reshape(mpix,peval.nx*peval.ny,peval.ncomp);
mwc=m.*wc; % pixels indicates locations and value of hte local maxima
mwcn=normcMax(mwc); % normalised to the maximum of each column. 