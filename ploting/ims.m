function h = ims(matrix, cmap, newfigure, noticks)
% function h = ims(matrix, cmap, newfigure, noticks)
% matrix - data to image
% cmap - colormap
% newfigure - creates new figure (default)
if nargin < 3 newfigure = 0; end
if nargin < 4 noticks = 0; end
    
matrixs = squeeze(matrix);
nd = ndims(matrixs);

if newfigure == 1; figure; end

switch nd
    case 2
        sd = sort(size(matrixs));
        if sd(1) == 1
            h = plot(matrixs);
        else            
            h = imagesc (matrixs);
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
            h(ii) = imagesc (matrixs(:,:,ii));
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