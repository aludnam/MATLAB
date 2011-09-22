function f = inverselinearfunction(x,p)
% f = inverselinearfunction(x,p)
%
% f = p(1)*1/x+p(2);
% p(1) = gradient
% p(2) = offset

f = p(1)*1./x+p(2);



