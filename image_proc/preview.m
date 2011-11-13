function out = preview(pathname,filename)
% preview(pathname,filename)
% reads one file (filename) and displays it. 
% (printmetadata.m prints some information from the metadata file)
if ~exist('pathname', 'var')
    pathname = '.';
end
if ~exist('filename','var')
    filename = 'img_000000001__000.tif';
    if ~exist(filename, 'file')
        filename='img_000000100_Sequence_000.tif';
    end
end
out = readim([pathname '/' filename]);
dipshow(out);