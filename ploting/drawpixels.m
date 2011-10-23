function drawpixels(sizeim,style, offset)
% drawpixels(sizeim,style, offset)
% draws pixels lines on the image (at integer + OFFSET positions)
% default style: ':'
% default offset: 0
if ~exist('style','var')
    style = ':';
end
if ~exist('offset','var')
    offset = [0 0];
end
xlim([0,sizeim(1)]+offset(1));
ylim([0,sizeim(2)]+offset(2));
for ii=0:sizeim(1)
    vline2(offset(1)+ii,style)
end
for ii=0:sizeim(2)
    hline(offset(2)+ii,style)
end
set (gca, 'DataAspectRatio',[1 1 1]);