function [x_n,y_n]=generatePointsOnLineX1X2(k,N,distrib,x1,x2,minDist)
% Generates N points disributed on a line of gradient k and between x-coordinates x1 and x2.
% [x_n,y_n]=generatePointsOnLine(k,N,distrib,X1,X2,minDist)
% x,y - coordinates
% k - gradient =tan(alpha)=opposite/adjecent
% N - number of points
% distrib -     'rand' uniform random distribution (defaults)
%               'regular' regular points
%               'randMinDist' readnom uniform but minimum distance larger
%               then mindist. Call:
%               [x,y]=generatePointsOnLine(k,N,distrib,mindist)
% Example:      
%               [x,y]=generatePointsOnLineX1X2(k,200,'regular',1,11);
%               scatter(x,y)
%               figure; 
%               [x_n,y_n]=generatePointsOnLineX1X2(.2,50,'randmindist',1,100,.8)
%               scatter(x_n,y_n)
if ~exist('minDist','var')
    minDist = 0;
end
[x,y]=generatePointsOnLine(k,N,distrib,minDist/(x2-x1));

x_n=x1+x*(x2-x1)*sqrt(1+k^2);
y_n=k*x_n;