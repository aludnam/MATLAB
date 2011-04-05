function out = normalize(in,p)
% out = normalize(in,p)
% normalizes 1d vector or 2d image or each slice of 3D image in the p-norm
% (sum(in.^p) = 1; (p=1 by default)
if nargin < 2 
    p=1;
end
nd=ndims(in);
switch nd
    case 1 %1D
        out = in/sum(in.^p);
    case 2 %2D
        out=bsxfun(@rdivide, in,sum(in.^p,1));
    case 3 %3D
        si = size(in);
        inr = reshape(in,si(1)*si(2),si(3));
        outr = bsxfun(@rdivide, inr, sum(inr.^p,1));
        out = reshape(outr, si);
end

