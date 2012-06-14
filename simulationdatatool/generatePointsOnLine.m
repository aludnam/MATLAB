function [x,y]=generatePointsOnLine(k,N,distrib,mindist)
% Generates N points disributed on a line of gradient k and length 1.
% [x,y]=generatePointsOnLine(k,N,distrib)
% x,y - coordinates
% k - gradient =tan(alpha)=opposite/adjecent
% N - number of points
% distrib -     'rand' uniform random distribution (defaults)
%               'regular' regular points
%               'randMinDist' readnom uniform but minimum distance larger
%               then mindist. Call: [x,y]=generatePointsOnLine(k,N,distrib,mindist)
% Example:      
%               [x,y]=generatePointsOnLine(0.5,100,'regular');
%               scatter(x,y)
% To make the line between x1 and x2 with gradient k and linear density r =
% N_n/(x2-x1)
% N_n=(x2-x1)*r;
% generate 
%               [x,y]=generatePointsOnLine(k,N_n,'regular');
% and 
%               x_n=x1+x*(x2-x1)*sqrt(1+k^2);
%               y_n=k*x_n;

if ~exist('distrib','var')
    distrib='rand';
end

q=1/sqrt(1+k^2);
switch lower(distrib)
    case 'rand'
        x = rand(N,1)*q;
    case 'regular'
        x = linspace(0,1,N)'*q;
    case 'randmindist'
        x = randMinDist(N,mindist); 
end
y = k*x;
end

function x=randMinDist(N,mindist)
ii=1;jj=1;
x=inf(N,1);
while ii<=N
    x(ii,1)=rand;      
    if min(abs(x(1:ii-1)-x(ii)))>mindist
        ii=ii+1;
    end   
    if ii==1 % x(1:ii-1)=Empty marix...
        ii=ii+1;
    end
    jj=jj+1;
    if or(jj>100000, N*mindist>1)
        error('Make the distance "mindist" smaller')
    end
end

end
