function h = ims(matrix, cmap, newfigure, noticks,clim)
% function h = ims(matrix, cmap, newfigure, noticks)
% matrix - data to image
% cmap - colormap
% newfigure - creates new figure (default)
% clim: [low high] is the intensity limits of the image. clim=[] fulls the colormap (default)
if nargin < 3; newfigure = 1; end
if nargin < 4; noticks = 0; end
if nargin < 5; clim = []; end
    
matrixs = squeeze(matrix);
nd = ndims(matrixs);

if newfigure == 1; figure; end

switch nd
    case 2
        sd = sort(size(matrixs));
        if sd(1) == 1
            h = plot(matrixs);
        else
            if ~isempty(clim)
                h = imagesc (matrixs,clim);
            else
                h = imagesc (matrixs);
            end
            set (gca, 'DataAspectRatio',[1 1 1]);
            if exist ('cmap', 'var')
                if ~isempty(cmap)
                    colormap(cmap);
                end
            end
        end
    case 3
        sm3 = size(matrixs,3);
        for ii=1:sm3
            if ii>1; figure; end
            if isempty(clim)
                h(ii) = imagesc (matrixs(:,:,ii));
            else
                h(ii) = imagesc (matrixs(:,:,ii),clim);
            end
            set (gca, 'DataAspectRatio',[1 1 1]);
            if exist ('cmap', 'var')
                if ~isempty(cmap)
                    colormap(cmap);
                end
            end
        end
end
if noticks
    set(gca,'xtick',[],'ytick',[])
end