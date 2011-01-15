function f = poissoncontunuos(x,lambda)
% function f = poissoncontunuos(x,lambda)
f = lambda.^x*exp(-lambda)./gamma(x+1);


