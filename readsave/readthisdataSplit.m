function [out, datadir] = readthisdataSplit(pathname, filename, range,nbatch)
% [out, datadir] = readthisdataSplit(pathname, filename, range,nbatch)
% Reads data. Splits reading so that it reades nbatch number of images at a time. It is faster (bug in dip_image?)
% Default for readimg data from Rainer and nbatch = 100;
% Example: [out, datadir] = readthisdataSplit('QD565_1',[],[0 100], 10);

if ~exist('nbatch', 'var')
    nbatch = 50;
end

nimages = range(2)-range(1)+1;
nsteps = ceil(nimages / nbatch);
sizeim = size(preview(pathname));

out = newim(sizeim(1),sizeim(2), nimages);

for ii=0:nsteps-1
    rangetmp(1) = range(1)+ii*nbatch;
    r = rangetmp(1)+nbatch-1;
    rangetmp(2) = min(r,range(2));    
    
    [outtmp, datadir] = readthisdata(pathname, filename, rangetmp);
    out(:,:,rangetmp(1):rangetmp(2)) = outtmp;
end