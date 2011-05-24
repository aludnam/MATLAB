function y=makeGauss(x,mu,sig)
% y=makeGauss(x,mu, sigma)
% Makes 1D gaussian centered at mu and with .
y = 1/sqrt(2*pi*sig^2)*exp(-(x-mu).^2/(2*sig^2));