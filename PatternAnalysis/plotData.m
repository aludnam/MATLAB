function plotData (data_in,box,colordata, sizedata)
% PLOTDATA plots 
% plotData (data_in,box,colordata, sizedata)
% data_in - N-by-2 matrix, each row is x-y coordinates of datapoint
% box - [xlim1 xlim2 ylim1 ylim2] - boundry of data
% colordata - color of ploted datapoints
% sizedata - size of ploted datapoints

if nargin<3 
    colordata = 'b';
end

if nargin<4 
    sizedata = 1;
end

if ~isempty(data_in)
    scatter(data_in(:,1), data_in(:,2),sizedata,colordata);
        
    if exist('box', 'var')
        if ~isempty(box)
            xlim([box(1), box(2)]); ylim([box(3) box(4)]);
        end
    end
end
set (gca, 'DataAspectRatio',[1 1 1]);
