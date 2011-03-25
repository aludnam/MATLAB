function [c1,c2] = cog(M)
% [c1,c2] = cog(M)
% Finds the center of gravity of the image M.

[rc,cc] = ndgrid(1:size(M,1),1:size(M,2));
Mt = sum(M(:));
c1 = sum(M(:) .* rc(:)) / Mt;
c2 = sum(M(:) .* cc(:)) / Mt;