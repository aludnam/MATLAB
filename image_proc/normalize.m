function out = normalize(in)
% out = normalize(in)
%normalizes 2d image
out = in/sum(in(:));
