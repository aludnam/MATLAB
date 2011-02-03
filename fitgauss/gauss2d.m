function G = gauss2d(sizevec, m, s, a)
% G = gauss2d(sizevec, m, s, a)
[X, Y] = meshgrid (1:sizevec(2), 1:sizevec(1));
G = a*exp(-((X-m(1)).^2+(Y-m(2)).^2)./(2*s^2));
