function d = ddivGauss(p, M)
% d = ddivGauss(p, M)
% D - divergence between M and gauss generated from p

G = gauss2d(size(M)', [p(1) p(2)], p(3),1)+1e-9;
% d = sum(sum((G - M).^2));
ddiv = ddivergence(G, M);
d = sum(ddiv(:));