function yfit=funkce(x,p)

omega = 2*pi/10;

yfit=p(1) + +p(2)*cos(omega*x).*exp(-p(3)*x.^2);