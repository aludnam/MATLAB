function f = poissoncontinuos(x,p)
% function f = poissoncontunuos(x,p)
% lambda = p(1)
% amplitude = p(2)
f = p(2)*p(1).^x*exp(-p(1))./gamma(x+1);
%  f = p(1).^x*exp(-p(1))./gamma(x+1);
%  f=f/sum(f);



