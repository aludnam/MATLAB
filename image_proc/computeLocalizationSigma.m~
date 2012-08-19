function sig=computeLocalizationSigma(lambda,na,pixelsize,background,nPhot)
% p=computeLocalizationSigma(lambda,na,pixelsize,background,nPhot)
% Computes std (sqrt(var)) of the localisation precision. 
% lambda: emisson wavelentgth
% na: numerical apperture
% pixelsize: size of the pixel
% background: background in one frame
% nPhot: number of detected photons
%
% See:  
% R. Thompson, "Precise Nanometer Localization Analysis for Individual Fluorescent Probes," Biophysical Journal 82, 2775?2783 (2002).
%
% Example: sig=computeLocalizationSigma(520,1.2,80,100,500)

ratio = lambda/na; 
sigAiry = sqrt(2)/(2*pi)*ratio; % Gaussian approximation of hte airy function
sig2 = (sigAiry^2+pixelsize^2/12)/nPhot + 8*pi*sigAiry^2*background^2/(pixelsize^2*nPhot^2); 
sig=sqrt(sig2); 
