function [xvec, yvec]=generatepointsregular(npoints,diamcirc,origin,phaseoffset)
% [xvec, yvec]=generatepointsregular(npoints,diamcirc,origin,pahseoffset)
% generates points qually distributed around the circle with 2r=diamcirc
% with center at origin
if ~exist('phaseoffset','var') %to make oblique configuration 
    phaseoffset = 0; 
end
angstep=2*pi/npoints;
angvectmp=(0:angstep:2*pi)+phaseoffset;
angvec=angvectmp(1:end-1); %last opint = first point
xvec=(origin(1)+0.5*diamcirc*cos(angvec))';
yvec=(origin(2)+0.5*diamcirc*sin(angvec))';
