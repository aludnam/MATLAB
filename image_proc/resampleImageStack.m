function out=resampleImageStack(in, zoom, interpolation_method)
% Uses resample.m function on the stack on images. Resamples each of the
% image in the scatck by zoom value.
% out=resampleImageStack(in, zoom, interpolation_method)
% default interpolation_method = '4-cubic';
% out(:,:,ii)=resample(in(:,:,ii),zoom, interpolation_method)
if ~exist('interpolation_method','var')
    interpolation_method = '4-cubic';    
    fprintf('Default interpolation method: %s\n',interpolation_method)
else
    fprintf('Interpolation method used: %s\n',interpolation_method)
end
nz=size(in,3); 

out=zeros(zoom*size(in,1),zoom*size(in,2), size(in,3)); 
for ii=1:nz
    out(:,:,ii)=resample(in(:,:,ii), zoom,0,interpolation_method);
end
