function out=subsampleconvolve(in, factor)
% out=subsampleconvolve(in, factor)

in = dip_image(in);
filter = ones(factor)/factor^2;
innew = convolve(in, filter);
out = resample(innew,1/factor);
