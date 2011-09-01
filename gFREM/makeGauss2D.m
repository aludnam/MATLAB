function y=makeGauss2D(x,mu,sig)
% y=makeGauss2D(x,mu, sigma)
% Makes 2D gaussian centered at [mu,0] and with variance sig^2 (diagonal and isotropic covariance matrix).

xx=x(:,:,1);
yy=x(:,:,2);

y = 1/(2*pi*sig)*exp(-((xx-mu).^2+yy.^2)/(2*sig^2));