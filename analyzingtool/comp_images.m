function comp_images(icapix, dpixc, d, rs, saveon, name, path)
if nargin < 7
    path = []; %saves to current directory
end
if nargin < 6
    saveon = 0;
end
nx = size(icapix,1)/rs;
ny = size(icapix,2)/rs;
[xm, ym, xfit,yfit] = companal(icapix, rs);
imstiled(imresize(icapix,1/rs),[],0);
hold off
if saveon
    SaveImageFULL([path 'components_' name], 'p')
end
sizedata = 40;

figure;
plotData([xfit; yfit]',[0 nx 0 ny],'xb',sizedata)
hold on;
plotData(d,[0 nx 0 ny],'r',sizedata)
hold off;
if saveon
    SaveImageFULL([path 'position_' name], 'p')
end

ims(imresize(sum(dpixc,3),1/rs),'gray');
hold on
plotData([xfit; yfit]',[0 nx 0 ny],'xb',sizedata);
hold off
if saveon
    SaveImageFULL([path 'comparioson_' name], 'p')
end