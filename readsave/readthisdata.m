function [out, datadir] = readthisdata(pathname, filename, range)
% [out, datadir] = readthisdata(pathname, filename, range)
% Reads data. Default for readimg data from Rainer.
% Example: [out, datadir] = readthisdata('QD565_1',[],[0 100]);

if ~exist('pathname', 'var')
    pathname = '.';
end
if or(~exist('filename','var'), isempty(filename))
        filename = 'img_000000*__000.tif';
end
if ~exist('range','var')
    range = [];
end

out = readtimeseries([pathname '/' filename],'',range);
dirnow = cd;
datadir = [dirnow '/' pathname];