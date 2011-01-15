function out = normalizeslices(in)
% out = normalizeslices(in)
%normalizes slices in 3d image

s=size(in);
inres = reshape(in, s(1)*s(2),s(3));
normvec=sum(inres,1);
normmat=repmat(normvec,s(1)*s(2),1);
outres=inres./normmat;
out = reshape(outres,s(1),s(2),s(3));
