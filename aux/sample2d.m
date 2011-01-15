function [x,y]=sample2d(image,n)
% [x,y]=sample2d(image,n);
% samples pixels from the distribution of the image
if ~exist('n', 'var')
    n=1;
end
i2=image-min(image(:)); %non negative but minimum 0;
i2n=normalize(i2);

xmarg =sum(i2n,2);
xmarg = xmarg/sum(xmarg(:));
xcdf=cumsum(xmarg);

ymarg =sum(i2n,1);
ymarg = ymarg/sum(ymarg(:));
ycdf=cumsum(ymarg);

rx=rand(1,n);
xxcdf=1:length(xcdf);
x=interp1(xcdf,xxcdf,rx);
ry=rand(1,n);
xycdf=1:length(ycdf);
y=interp1(ycdf,xycdf,ry);

% 
% d=diff(abs(xcdf-rand)); 
% [v,x]=max(abs(d));
% d=diff(abs(ycdf-rand)); 
% [v,y]=max(abs(d));


% [v, x]=min(abs(xcdf-rand));
% [v, y]=min(abs(ycdf-rand));