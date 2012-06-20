function p=computeAirySigmaFWHM(lambda, na)
% p=computeAirySigmaFWHM(lambda, na)
% Computed radius of the airy disk (p.airy), std of the Gaussian
% approximation (p.sigma), adn FWHM of hte Gaussian approximation (p.FWHM)
% from the wavelength (lambda) and numerical apperture (na); 
%
% Example: p=computeAirySigmaFWHM(520, 1.2)

p.ratio = lambda/na; 
p.airy = .61*p.ratio; % radius of the airy disk
p.sigma = sqrt(2)/(2*pi)*p.ratio; % Gaussian approximation of hte airy function
p.FWHM = 2*sqrt(2*log(2))*p.sigma; % full width in half maximum of the gauss approximation