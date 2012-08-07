function [onesource,twosources,edgedist,nlocmax]=classifyW_localmin(w,nx,ny,psf,threshold,displayimage)
% [onesource,twosources,edgedist,nlocmax]=classifyW_localmin(w,nx,ny,psf,threshold,displayimage)
% Classify the the estimated sources w (nx*ny X ncomp matrix). 
%
% Input: 
% w - (nx*ny X ncomp) matrix of NMF estimated sources
% nx,ny -> dimensions of the 2D image of the sources
% threshold -> fraction of the global maximum of each w_k above which the
% local miminum must be to be considered for evaluation (default: threshold=.5)
% psf -> image of the point spread function OR struct containing field psf.NA,psf.lambda,psf.pixelisze
% displayimage -> display the image of classified sources. 
%
% Output: 
% onesource -> binary vector indicating 1 source
% twosources -> binary vector indicating 2 sources
% edgedist  -> distance of the local maximum from the edge
% nlocmax -> numbe of local maxima in each 
%
% Example: 
% [onesource,twosources,edgedist,nlocmax]=classifyW_localmin(res.w(:,1:end-1),peval.nx,peval.ny,p,0.5,1)

if isstruct(psf)    
    psf = PSFGEN('lambda', psf.lambda, 'na', psf.NA, 'pixelsize', psf.pixelsize, 'sizevec', [nx ny], 'method', 'airy', 'nphot',1, 'verbose',0);
    fprintf('Generating PSF...\n')
end

if ~exist('threshold','var')
    threshold = 0.5;
end

if ~exist('displayimage','var')
    displayimage = 0;
end

[mwcn,wpixc,globmax,edgedist]=localminW(w,nx,ny,psf);
nlocmax = sum(mwcn>threshold); % number of local maxima, thresholed by 50% of the global maximum
onesource=nlocmax==1;
twosources=nlocmax==2;

if displayimage
    framecolvec=(edgedist'<2)+(2*onesource);
    framecolvec(framecolvec==1)=0;
    framecolvec=framecolvec+twosources;
    figure; imstiled(reshape(w,nx,ny,size(w,2)),[],'gray',[1:size(w,2)],[],[],framecolvec)
end
