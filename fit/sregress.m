function [p,sp] = sregress(x,y,sy,n)
%
%  [p,sp] = regress(x,y,sy,n) fits a polynomial of degree n to data
%   points (x,y). sy contains 1 sigma errors in y, and sp returns the
%   1 sigma errors in the polynomial coefficients p. The best polynomial
%   is  y = p(1) x^n + p(2) x^(n-1) + ... + p(n-1)
%
x = x(:); y = y(:); sy = sy(:);
lx = length(x);
A = zeros(lx,n+1);
% Construct Vandermonde matrix
for j=1:n+1
        A(:,j) = x.^(n+1-j);
end
AWA = A'*diag(1 ./sy.^2)*A;
AWy = A'*diag(1 ./sy.^2)*y;
AWAi = inv(AWA);
p = AWA\AWy;
sp = sqrt(diag(AWAi));
