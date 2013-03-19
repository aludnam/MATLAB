function G = gauss2d(sizevec, m, s, a)
% G = gauss2d(sizevec, m, s, a)
% Gaussian function of size given by sizevec(1)Xsizevec(2) with position
% vecotr m, standard deviation s and intensity a. 
[X, Y] = meshgrid (1:sizevec(2), 1:sizevec(1));
G = a*exp(-((X-m(1)).^2+(Y-m(2)).^2)./(2*s^2));
