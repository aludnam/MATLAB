function n = normcSum(m)
%NORMC Normalize columns of a matrix.
%
%  Syntax
%
%    normcSum(M)
%
%  Description
%
%    NORMCSUM(M) normalizes the columns of M to a sum of 1.
%
%  Examples
%    
%    m = [1 2; 3 4]
%    n = normc(m)
%
%  See also NORMC

[mr,mc] = size(m);

n = m./repmat(sum(m),mr,1);