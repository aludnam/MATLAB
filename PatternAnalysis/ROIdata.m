function [dataXY_out] = ROIdata (dataXY_in, x_lim1, x_lim2, y_lim1, y_lim2, plotOrig)
% [dataXY_out] = ROIdata (dataXY_in, x_lim1, x_lim2, y_lim1, y_lim2,
% plotOrig)
% reduction of the data size
% plotOrig - plot original data with selected boundary

if nargin<6
    plotOrig = 0; %not plotting by default
end
x_in = dataXY_in (:,1);
y_in = dataXY_in (:,2);
x_out = [];
y_out = [];

jj = 1;

sizeData = length(x_in);
for ii = 1:sizeData
    if and(x_in(ii) > x_lim1, x_in(ii) <= x_lim2)
        if and(y_in(ii) > y_lim1, y_in(ii) <= y_lim2)
            x_out(jj) = x_in(ii);
            y_out(jj) = y_in(ii);
            jj = jj+1;
        end
    end
end
            
if plotOrig
    figure
    scatter(x_in, y_in, 1);
    %draw box:
    line([x_lim1 x_lim2],[y_lim1 y_lim1], 'Color','r', 'LineWidth', 3);
    line([x_lim2 x_lim2],[y_lim1 y_lim2], 'Color','r', 'LineWidth', 3);
    line([x_lim1 x_lim2],[y_lim2 y_lim2], 'Color','r', 'LineWidth', 3);
    line([x_lim1 x_lim1],[y_lim1 y_lim2], 'Color','r', 'LineWidth', 3);
    set (gca, 'DataAspectRatio',[1 1 1]);
end

dataXY_out = [];
if ~isempty(x_out)
    dataXY_out (:,1) = x_out';
    dataXY_out (:,2) = y_out';
end