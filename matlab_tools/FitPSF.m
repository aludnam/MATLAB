function [FWHM]=FitPSF(img)
%[fittedSum,paramsSum,idiv,myfunct] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))+c(5)',[(max(img)-min(img)) 1 1 8 min(img)],img,1000);
%FWHM=sqrt(paramsSum(4))* 2*sqrt(log(2));  % to get FWHM
[paramsSum,mse,fittedSum] = FitDataNDFast([min(img) 5 (max(img)-min(img)) 0 0],img,2,1000,'mse');
FWHM=sqrt(paramsSum(2))* 2*sqrt(log(2));  % to get FWHM

img-fittedSum

ypos=floor(size(fittedSum,2)/2)+round(paramsSum(4));
%ypos=floor(size(fittedSum,2)/2)+round(paramsSum(3));

figure;
plot(fittedSum(:,ypos),'b');hold on; plot(img(:,ypos),'r');
%dipshow(11,img-fittedSum); drawnow;
