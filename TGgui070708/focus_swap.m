% Handles experimental focus swapping

function focus_swap(hObject, eventdata, handles)
    x_size     = str2double(get(handles.x_size_edit,'String'));
    y_size     = str2double(get(handles.y_size_edit,'String'));
    
    if x_size < y_size
        positionf = [.18,.065,.87,.87];
    else
        positionf = [.32,.3,.6,.65];
    end

    focused = getappdata(0,'focused');

    if focused == hObject %Return when clicking an object already in focus
        return;
    end

    current = getappdata(hObject,'pos');

    set(focused, 'Position', current);      %Swap the axes positions
    axis image
    setappdata(focused,'pos',current);

    set(hObject, 'Position', positionf);
    axis image
    setappdata(0,'focused', hObject);


    if get(handles.show_colorbar,'Value')
        pos = get(hObject, 'Position');
        colorbar('Position',[.947, pos(2), .015, pos(4)],'DrawMode','fast')
    end

    axes(focused)
    colorbar('delete');
    drawnow