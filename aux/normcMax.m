function n = normcMax(m)
%NORMC Normalize columns of a matrix such taht the maximum is 1.
%
%  Syntax
%
%    normcSum(M)
%
%  Description
%
%    NORMCSUM(M) normalizes the columns of M to a max of 1.
%
%  Examples
%    
%    m = [1 2; 3 4]
%    n = normc(m)
%
%  See also NORMC

[mr,mc] = size(m);

n = m./repmat(max(m,[],1),mr,1);