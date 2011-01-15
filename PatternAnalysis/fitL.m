function y = fitL(x,p)
% function y = fitL(x,p)
% fits L function for Thomas Process
% p(1) - kappa, p(2) - sigma
y = sqrt(x.^2 + 1/(p(1)*pi) * (1-exp(-(x.^2/(4*p(2)^2))))) - x;