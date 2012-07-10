function out = binsumcol(in, numinbin)
% Rebins image by summing numinbin bins together along columns of the
% 2D data.
% out = binsumcol(in, numinbin)

di=0;
if strcmp(class(in),'dip_image')
    di=1;
    in=double(in);
end
tmp = in;

for ii=2:numinbin
    tmp = tmp + circshift(in,[-1*(ii-1) 0]);
end

nb = floor(size(in,1)/numinbin);
ix=numinbin*(1:nb)-(numinbin-1);
out = tmp(ix,:);
if di
    out=dip_image(out);
end