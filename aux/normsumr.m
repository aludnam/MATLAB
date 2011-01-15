function out = normsumr(in)
% out = normsumr(in)
% nomelizes sum of rows to 1

[m,n]=size(in);
sumin = sum(in,2);
out = in./repmat(sumin,1,m); %normalization of each component