function [xout,yout]=rejectsample(f, n)
% [xout,yout]=rejectsample(f, n)
% rejection sampling from the 2D discrete distribution f 
% n: number of samples
% the proposal distriburion is a uniform distribution on the range of f
sizef=size(f);
maxf=max(f(:));
nsofar=0;
while nsofar<n
    xrand=sizef(1)*rand;
    yrand=sizef(2)*rand;
    u=maxf*rand;
    if u<f(ceil(xrand), ceil(yrand))
        nsofar=nsofar+1;
        x(nsofar)=xrand;
        y(nsofar)=yrand;
    end    
end

%to the middle of the pixels...
xout=x-0.5;
yout=y-0.5;
    
