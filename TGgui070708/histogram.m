%This function plots a histogram, given a -grayscale- file input
%NOTE: *** Possible errors with 32 bit, since this program doesnt read the alpha channel ***
%
%Last Modified 8/15/07

function [maximum] = histogram(handles, filename, bit_depth, index)
    if nargin == 1
        cur = pwd; %Find the current directory
        if isappdata(0, 'Last_Directory')
            last_dir = getappdata(0, 'Last_Directory');
            if last_dir ~= 0
                if exist(last_dir,'dir')
                    cd(last_dir);
                end
            end
        end

        [fname,pname] = uigetfile({'*.tif;*.tiff;*.bmp;*.jpg;*.jpeg;*.png;*.gif','Image Files (*.tif, *.bmp, *.jpg, *.png, *.gif)'},'Select a grayscale image file:');

        if pname ~=0
            if exist(pname,'dir')
                setappdata(0, 'Last_Directory', pname);
            end
        end
        cd(cur); %Reset active directory

        if fname == 0
            return
        else
            filename = [pname,fname];
        end

        file_info = get_file_info2(filename, handles, 1);

        if strcmp(file_info.type,'multipage')
            x = inputdlg(sprintf('Multi-page image file detected.\n\nEnter index of image:'),'',1,{'1'});
            if ~isempty(x)
                image = imread(filename,str2double(x));
            else
                return
            end
        else
            image = imread(filename);
        end

        bit_depth = file_info.bitdepth;

        n_bins = ceil((2^bit_depth)/10);%round(2*(x_size*y_size)^(1/3));
        bins = linspace(0, 2^bit_depth, n_bins);

        figure('Name','Histogram')
        hist(image(:), bins)            %Create histogram of image
        h_hist = findobj(gca,'Type','patch');
        set(h_hist,'EdgeColor','k');
                                        %Set title to filename with no TeX support
        if ~exist('x','var')
            title(fname,'Interpreter','none');
        else
            title(sprintf([fname,'\nIndex:',char(x)]),'Interpreter','none');
        end

        xlabel(['Image Intensity (grayscale, 0-',num2str(2^bit_depth-1),')']);     %Set appropriate labels
        ylabel('Counts');
    elseif nargin > 1
        if nargin == 3
            image = imread(filename);
        elseif nargin == 4
            image = imread(filename,index);
        else
            error('Error: Invalid number of input arguments.');
        end

        n_bins = ceil((2^bit_depth)/10);
        bins = linspace(0,2^bit_depth,n_bins);

        h = hist(image(:),bins);
        figure('NumberTitle','off','Name','Histogram: Auto-selected Background Level (red)');
        [n xout] = hist(image(:),bins);

        [c,i] = max(h);

        %Ignore peaks at 0, this is likley the wrong peak
        if i == 1  %1=>0 in this case
            [c,i] = max(h(2:size(h,2)));
            i=i+1;
        end

        nn = n(i);
        xn = xout(i);

        hold on
        n(i) = 0;
        bar(xout, n, 'BarWidth',1);
        bar(xn, nn,'FaceColor','r','BarWidth',10);
        a = get(gcf,'CurrentAxes');
        set(a,'XLim',[xn-500 xn+1700],'YLim',[0 nn+.05*nn]);

        [path,name] = fileparts(filename);
        if exist('index','var')
            title(['*\',name,' -- Index:',num2str(index)],'Interpreter','none');
        else
            title(['*\',name],'Interpreter','none');
        end

        hold off
        drawnow

        maximum = bins(i);
    end