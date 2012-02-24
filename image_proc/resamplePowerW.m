function wout=resamplePowerW(w,nx,ny,resampleFactor,power,makeImage)
% Resamples (by resampleFactor) resutlts W of NMF (res.w), makes their power (by power) and normalizes
% them to sum to 1;
% makeImage = 1 results in wout being stack of images (nx X ny X ncomp)
% rather then a nx*ny X ncomp matrix (default)
%
% wout=resamplePowerW(w,nx,ny,resampleFactor,power,makeImage)
%
% example: 
% wout=resamplePowerW(res.w, peval.nx, peval.ny, 4,2,1);
% figure; imstiled(wout)

if ~exist('makeImage','var')
    makeImage = 0; 
end
ncomp=size(w,2); 
wr=resampleImageStack(reshape(w, nx, ny, ncomp),resampleFactor).^power;
wout=normcSum(reshape(wr,resampleFactor^2*nx*ny,ncomp)); % normalized to sum to 1
if makeImage
    wout = reshape(wout, resampleFactor*nx, resampleFactor*ny, ncomp); 
end