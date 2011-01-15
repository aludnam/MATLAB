% Executes on pause button press, is its own file so that it can
% be used globally throughout
function pause_gui(hObject, eventdata)
    if get(hObject,'UserData')
        set(hObject,'String','Resume','UserData',0);
        uiwait;
    else
        set(hObject,'String','Pause','UserData',1);
        uiresume;
    end