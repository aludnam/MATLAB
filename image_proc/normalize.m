function out = normalize(in)
% out = normalize(in)
%normalizes 1d vector or 2d image or each slice of 3D image
nd=ndims(in);
switch nd
    case 1 %1D
        out = in/sum(in);
    case 2 %2D
        out=bsxfun(@rdivide, in,sum(in,1));
    case 3 %3D
        si = size(in);
        inr = reshape(in,si(1)*si(2),si(3));
        outr = bsxfun(@rdivide, inr, sum(inr,1));
        out = reshape(outr, si);
end

