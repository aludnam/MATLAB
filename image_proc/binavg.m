function out = binavg(in, numim2bin)
% out = binavg(in, numim2bin)
% rebins image by averaging nnumim2bin bins together

di=0;
if strcmp(class(in),'dip_image')
    di=1;
    in=double(in);
end
tmp = in;

for ii=1:numim2bin
    tmp = tmp + circshift(in,[0 0 -1*ii]);
end

nb = floor(size(in,3)/numin2bin);
ix=numim2bin*(1:nb)-(numin2bin-1);
out = tmp(ix)/numin2bin;
if di
    out=dip_image(out);
end