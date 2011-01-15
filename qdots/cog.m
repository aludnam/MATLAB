function [c1,c2] = cog(M)
% center of gravity

[rc,cc] = ndgrid(1:size(M,1),1:size(M,2));
Mt = sum(M(:));
c1 = sum(M(:) .* rc(:)) / Mt;
c2 = sum(M(:) .* cc(:)) / Mt;