function d = difference(p, M)

G = gauss2d(size(M)', [p(1) p(2)], p(3),1)+1e-9;
% d = sum(sum((G - M).^2));
ddiv = ddivergence(G, M);
d = sum(ddiv(:));