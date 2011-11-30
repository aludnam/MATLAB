function out = binsum(in, numim2bin)
% Rebins image by averaging nnumim2bin bins together along z-direction
% out = binavg(in, numim2bin)


if strcmp(class(in),'dip_image')
    di=1;
    in=double(in);
end
tmp = in;

for ii=1:numim2bin-1
    tmp = tmp + circshift(in,[0 0 -1*ii]);
end

nb = floor(size(in,3)/numim2bin);
ix=numim2bin*(1:nb)-(numim2bin-1);
out = tmp(:,:,ix);