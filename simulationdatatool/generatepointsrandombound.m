function [xvec, yvec]=generatepointsrandombound(npoints,diamcirc,origin)
% [xvec, yvec]=generatepointsrandombound(npoints,diamcirc,origin)
% generates points randomly distributed in the circle with 2r=diamcirc
% with center at origin

xvec=(origin(1)+diamcirc*(rand(1,npoints)-.5))';
yvec=(origin(2)+diamcirc*(rand(1,npoints)-.5))';
