function out=shiftmulti(in, shiftmatrix)
% out=shiftmulti(in, shiftmatrix)

sv = size(in);
out = zeros(sv);
for ii=1:sv(3)
    out(:,:,ii)=shift(in(:,:,ii), shiftmatrix(ii,:));
end