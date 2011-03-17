function [out, meanvalues]=makezeromean(in,dim)
% [out, meanvalues]=makezeromean(in,dim)
% Makes IN zero mean along the dimension DIM. Does not alter the variance.

meanvalues=mean(in,dim);
out=bsxfun(@minus,in,meanvalues);