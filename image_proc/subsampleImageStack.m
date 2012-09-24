function out=subsampleImageStack(in, subsample_factor)
% Uses subsample.m function on the stack on images. Subsample each of the
% image in the scatck by zoom value.
% out=subsampleImageStack(in, subsample_factor)
% default subsample_factor=2
% out(:,:,ii)=subsample(in(:,:,ii),subsample_factor)
if ~exist('subsample_factor','var')
    subsample_factor=2;
    fprintf('Default subsample factor: %d\n',subsample_factor)
end
   
nz=size(in,3); 

out=zeros(size(in,1)/subsample_factor,size(in,2)/subsample_factor, size(in,3)); 
in=double(in); 
for ii=1:nz
    out(:,:,ii)=double(subsample(squeeze(in(:,:,ii)), subsample_factor))';
end
