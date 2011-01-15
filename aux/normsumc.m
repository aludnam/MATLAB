function out = normsumc(in)
% out = normsumc(in)
% nomelizes sum of colum to 1

[m,n]=size(in);
sumin = sum(in,1);
out = in./repmat(sumin,n,1); %normalization of each component