% Combines image files when given a list, these can be multipage files or
% single page files.
%
% 07/17/07

function combine(path, path_out, filelist)
% tic
    if iscellstr(filelist)
        num = numel(filelist);
        filelist = sort(filelist);  %Sort by ascii name

        w = waitbarxmod(.01,'Combining...','CreateCancelBtn','delete(gcf)');    %Create waitbar
        set(w,'Name','Progress Bar');
        pause(.1);
        keepontop('Progress Bar');

        for loop=1:num
            info = imfinfo([path,filelist{loop}]);
            lnth = length(info);
            infile = [path,filelist{loop}];
%             current_file = filelist{loop}

            for idx=1:lnth
                image0 = imread(infile,idx);
                imwrite(image0,path_out,'Compression','none','WriteMode','append');
                waitbarxmod(loop/num,w);
                drawnow;
            end
        end
    else
        return
    end

% total_time_sec = toc

    if exist('w','var')
        waitbarxmod(1);
        drawnow;
        pause(.1);

        delete(w);
    end