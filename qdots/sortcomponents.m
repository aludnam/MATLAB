function [sx, isx, fsx] = sortcomponents(x,p)
% [sx, isx, fx] = sortcomponents(x, p)
% Sorts the inut values x (column vectores - like res.w) according to the
% x.^p score. 
% p: exponent for the score. (p-norm) default p=2; 
% sx: sorted compontnts;
% isx: index of sorting - x(:,isx) are sorted...
% fsx: sorted objective function used for sorting
% fx = mean (x,2);
% ts = testimportance(reshape(dpixc, peval.nx*peval.ny, peval.nt), res.w, res.h);

if ~exist('p','var')
    p=2;
end
fx= sum(x.^p,1);
[sx, isx] = sort(fx, 'descend');
fsx=fx(isx);
