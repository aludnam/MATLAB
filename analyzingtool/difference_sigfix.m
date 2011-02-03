function d = difference_sigfix(p, M, sigfix)

G = gauss2d(size(M), [p(1) p(2)], sigfix,1);
d = sum(sum((G - M).^2));
