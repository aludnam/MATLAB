function [xvec, yvec]=generatepointsrandombound(npoints,diamcirc,origin)
% [xvec, yvec]=generatepointsrandombound(npoints,diamcirc,origin)
% generates points randomly distributed in the circle with 2r=diamcirc
% with center at origin

xvec=(origin(1)+0.5*diamcirc*rand(1,npoints))';
yvec=(origin(2)+0.5*diamcirc*rand(1,npoints))';
