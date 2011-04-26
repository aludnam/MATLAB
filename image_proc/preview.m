function preview(pathname,filename)
% preview(pathname,filename)
% reads one file (filename) and displays it. 
if ~exist('pathname', 'var')
    pathname = '.';
end
if ~exist('filename','var')
    filename = 'img_000000001__000.tif';
end
out = readim([pathname '/' filename]);
dipshow(out);