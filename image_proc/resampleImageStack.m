function out=resampleImageStack(in, zoom)
% Uses resample.m function on the stack on images. Resamples each of the
% image in the scatck by zoom value.
% out=resampleImageStack(in, zoom)
% out(:,:,ii)=resample(in(:,:,ii),zoom)
nz=size(in,3); 
out=zeros(zoom*size(in,1),zoom*size(in,2), size(in,3)); 
for ii=1:nz
    out(:,:,ii)=resample(in(:,:,ii), zoom,0,'4-cubic');
end
