%get_file_info2.m
%
%Replacement for get_file_info.m:
%   Numerous efficiency fixes, smaller filesize, easier readability,
%   and structural/organizational changes to keep returned info together

function [file_info] = get_file_info2(file, handles, override)
    [file_info.path, file_info.name, file_info.ext] = fileparts(file);

    try
        image0    = imread(file,2);
        file_info.type = 'multipage';
        int_prec=0;
    catch
        image0    = imread(file);
        file_info.type = 'singlepage';
    end

    class_type    = class(image0);
    file_info.bitdepth = log2(double(intmax(class_type))-double(intmin(class_type)));
    file_info.width    = size(image0,2);
    file_info.height   = size(image0,1);

    if ~isnan(str2double(get(handles.prec_edit,'String')))
        int_prec = str2double(get(handles.prec_edit,'String'));
        file_info.prec = ['%',num2str(int_prec),'.',num2str(int_prec),'d'];
    else
        length_name = length(file_info.name);
        for x = length_name:-1:1
            if isnan(str2double(file_info.name(x)))
                int_prec  = length_name-x;
                file_info.prec = ['%',num2str(int_prec),'.',num2str(int_prec),'d'];
                break
            end
        end
    end

    file_info.part_name = file_info.name(1:length(file_info.name)-int_prec);

    if get(handles.custom_frame,'Value')
        file_info.start = str2double(get(handles.from_edit,'String'));
        file_info.stop  = str2double(get(handles.to_edit,'String'));
    elseif ~get(handles.custom_frame,'Value')
        if strcmp(file_info.type,'multipage')
            file_info.start = 1;

            if ~exist('override','var') %Allow function inputs to override autodetect of end frame (saves on time)
                file_imf = imfinfo(file);
                file_info.stop  = length(file_imf);
            end
        elseif strcmp(file_info.type,'singlepage')
            file_info.start = 1;

            if ~exist('override','var') %Allow function inputs to override autodetect of end frame (saves on time)
                index = 1;
                while exist([file_info.path,'\',file_info.part_name,num2str(index,file_info.prec),file_info.ext],'file')
                    index = index + 1;
                end
                file_info.stop = index - 1;
            end
        end
    end

    if strcmp(file_info.type,'multipage') %Override prec and partial name for multipage images
        file_info.part_name = file_info.name;
        file_info.prec  = '%0.0d';
    end