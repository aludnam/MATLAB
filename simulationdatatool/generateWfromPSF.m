function w=generateWfromPSF(sizeVec,x_vec,y_vec,psf)
% w=generateWfromPSF(sizeVec,x_vec,y_vec,psf)

npoints = length(x_vec);
center=sizeVec/2;
centerPsf=floor([size(psf,1),size(psf,2)]/2);
shift_y=x_vec-center(1); 
shift_x=y_vec-center(2);
w=zeros(prod(sizeVec),npoints);
if ndims(psf)<3
    % if pased only one 2D image of PSF then all points have the same PSF
    psf = repmat(psf,[1,1,npoints]);
end

for ii=1:npoints
    wpix=shift(squeeze(psf(:,:,ii)),[shift_x(ii),shift_y(ii)]);
    lVec=centerPsf-floor(center);    
    wpixS=wpix(lVec(1):lVec(1)+sizeVec(1)-1,lVec(2):lVec(2)+sizeVec(2)-1);
    w(:,ii)=reshape(wpixS,prod(sizeVec),1);
end