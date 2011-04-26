function out=catImageMeanVar(in)
% out=catImageSumVar(in)
%
% Concatenate image (in) with mean(in,[],3) and var(in,[],3) into one large
% image
if ~strcmp(class(in),'dip_image')
    in = dipimage(in);
end
sv = size(in,3);
dold=dipgetpref('DefaultGlobalStretch');
dipsetpref('DefaultGlobalStretch', 'OFF');
out = cat(4, in,repmat(mean(in,[],3),1,1,sv), repmat(var(in,[],3),1,1,sv))
dipsetpref('DefaultGlobalStretch',dold);