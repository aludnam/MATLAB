function [xvec, yvec]=generatepointsrandomboundCircle(npoints,diamcirc,origin, minsep)
% [xvec, yvec]=generatepointsrandomboundCircle(npoints,diamcirc,origin,
% minsep)
% generates points randomly distributed in the circle with 2r=diamcirc
% with center at origin. Minimum separation between points is minsep. 
if ~exist('minsep','var')
    minsep = 0;
end

n=1;
xvec = 0;
yvec = 0;
while n<=npoints
    pos=origin+diamcirc*(rand(1,2)-.5);
    if sum((pos-origin).^2) <= (diamcirc/2)^2
        if or(n==1, min(pdist([[xvec'; pos(1)],[yvec'; pos(2)]], 'euclidean'))>minsep)
            xvec(n)=pos(1);
            yvec(n)=pos(2);
            n=n+1;
        end
    end
end
