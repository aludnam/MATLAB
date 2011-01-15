function y=gauss_res2(p)
global grab;
global xpix;
global ypix;
global psf_w02;

yfit=p(3)*exp(-2*((xpix-p(1)).*(xpix-p(1))+(ypix-p(2)).*(ypix-p(2)))/psf_w02);
ydev=yfit-double(grab);
ydev2=ydev.*ydev;
y=sum(sum(ydev2));
%y=100*(p(2)-p(1)^2)^2+(1-p(1))^2;
end
