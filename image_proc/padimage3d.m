function out=padimage3d(in,padsize, padvalue)
% out=padimage3d(in,padsize)
% Pads 3d image (in) with padvalues (default: 0) into a new 3d image (out) of the size [padsize(1), padsize(2), size(in,3)]
if ~exist('padvalue','var')
    padvalue = 0;
end
sin = size(in);
out = zeros(padsize(1), padsize(2), sin(3));
for ii=1:sin(3)
    impad=padimage(in(:,:,ii),padsize);
    impad(impad==0)=padvalue;
    out(:,:,ii)=impad;
end