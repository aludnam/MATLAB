function [timeson, timesoff, fbin]=ontimestat(bm, threshold)
% [timeson, timesoff, fbin]=ontimestat(bm, thershold)
% Computes ON and OFF times of the time serires bm. ON (OFF) time is how long is
% the time-serires bm above (below) the threshold.
% bm=f';
% threshold=500;

fbin=(bm>threshold);
timeson=computetimes(fbin);
timesoff=computetimes(1-fbin);
end
function times=computetimes(fbin)
zv = zeros(1,size(fbin,2));
fbin2=[zv; fbin; zv];
pos=find(diff(fbin2)>0);
neg=find(diff(fbin2)<0);
times=neg-pos;
end