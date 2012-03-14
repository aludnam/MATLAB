function [L,nout]=bwlabelStack(bw,n)
% See BWLABEL help -> this is just made for a stack of images.
z=size(bw,3);
for ii=1:z;
    [L(:,:,ii),nout(ii)] = bwlabel(squeeze(bw(:,:,ii)), n);
end
    