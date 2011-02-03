function out = normalize(in)
% out = normalize(in)
%normalizes 1d vector or 2d image or each slice of 3D image
nd=ndims(in);
if nd<3 %1D or 2D
    out = in/sum(in(:));
else %3D
    si = size(in);
    inr = reshape(in,si(1)*si(2),si(3));
    outr = bsxfun(@rdivide, inr, sum(inr,1));
    out = reshape(outr, si);
end
    
