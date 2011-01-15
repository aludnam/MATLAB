%A function to handle gui errors

function rethrow_gui_err(err)
    if nargin == 0
        err = lasterror;
    end

    if isfield(err,'stack')
        msgbox([sprintf(['Line: ',num2str(err.stack(1).line),'\nFunction: ',err.stack(1).name,'\n\n']),err.message],'Error','error','replace')
    else
        msgbox('An error has occurred: See the command window for details of the last known error.','Error','error','replace');
    end

    rethrow(lasterror) %Rethrow the last error and return control to the keyboard