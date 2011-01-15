function f = gaussfunction(x,p)
% f = gaussfunction(x,p)
% p(1) = mean
% p(2) = std deviation (sqrt(var))
% p(3) = scale coefficient (if p(3)==1 then it is normalised)
normconst = 1/(sqrt(2*pi*p(2)^2));
f = p(3)*normconst*exp(-((x-p(1)).^2)/(2*p(2)^2));



