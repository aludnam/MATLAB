%   Copyright 2008, S. T. Hess Lab, University of Maine, Orono, ME
%   Do not distribute without permission
%
%   gui_cp.m
%   ---------------------------------------------------
%   Displays best with desktop resolution 1024x768 and higher
%
%   Note that only 1 function should be running at a time. Multiple
%   functions running (multiple analysis button presses) will cause errors 
%   since matlab does -not- support multithreading.

function varargout = gui_cp(varargin)
%------------------------------------------------------
%DO NOT EDIT - Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @gui_cp_OpeningFcn, ...
                       'gui_OutputFcn',  @gui_cp_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
%DO NOT EDIT - End initialization code - DO NOT EDIT
%------------------------------------------------------


% --- Executes just before gui is made visible.
function gui_cp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    if isappdata(0,'last_mbkg_save')   %Clear appdata that will be used,
        rmappdata(0,'last_mbkg_save'); %checks for existence first so that script doesnt error out
    end                                %before the GUI loads
    if isappdata(0,'last_einzel_save')
        rmappdata(0,'last_einzel_save');
    end
    if isappdata(0,'save')
        rmappdata(0,'save');
    end
    if isappdata(0,'temp_runtime')
        rmappdata(0,'temp_runtime');
    end
    if isappdata(0,'from_mat_file');
        rmappdata(0,'from_mat_file');
    end

    set(handles.roi_pre,'Color',[.941,.941,.941],'Box','on');

%Setup custom GUI tabs (low level)               |
%---------------------------------------         |
%                                                v
    axes(handles.localization_toggle);
   	text('String','Localization',...
         'Units','Normalized',...
         'Position',[.5,.5],... %Position in center of axes
         'HorizontalAlignment','center',...
         'VerticalAlignment','middle',...
         'Margin',0.001,...
         'FontSize',9,...
         'FontName', 'MS Sans Serif',...
         'HitTest','off'); %HitTest off so text doesnt block ButtonDownFcns of the axes
%                               'Backgroundcolor',[.753, .753, .753],...
%                               'ButtonDownFcn','gui_cp(''threshold_toggle_Callback'',gcbo,[],guidata(gcbo))');

    axes(handles.background_toggle);
    text('String','Background',...
         'Units','Normalized',...
         'Position',[.5,.5],...
         'HorizontalAlignment','center',...
         'VerticalAlignment','middle',...
         'Margin',0.001,...
         'FontSize',9,...
         'FontName', 'MS Sans Serif',...
         'HitTest','off');

    axes(handles.roi_toggle);
    text('String','ROI',...
         'Units','Normalized',...
         'Position',[.5,.5],...
         'HorizontalAlignment','center',...
         'VerticalAlignment','middle',...
         'Margin',0.001,...
         'FontSize',9,...
         'FontName', 'MS Sans Serif',...
         'HitTest','off');

    axes(handles.optics_toggle);
    text('String','Optics',...
         'Units','Normalized',...
         'Position',[.5,.5],...
         'HorizontalAlignment','center',...
         'VerticalAlignment','middle',...
         'Margin',0.001,...
         'FontSize',9,...
         'FontName', 'MS Sans Serif',...
         'HitTest','off');

    axes(handles.preview_toggle_tab);
    text('String','Preview',...
         'Units','Normalized',...
         'Position',[.5,.5],...
         'HorizontalAlignment','center',...
         'VerticalAlignment','middle',...
         'Margin',0.001,...
         'FontSize',9,...
         'FontName', 'MS Sans Serif',...
         'HitTest','off');
     
    axes(handles.render_toggle);
    text('String','Render',...
         'Units','Normalized',...
         'Position',[.5,.5],...
         'HorizontalAlignment','center',...
         'VerticalAlignment','middle',...
         'Margin',0.001,...
         'FontSize',9,...
         'FontName', 'MS Sans Serif',...
         'HitTest','off');
%---------------------------------------

    if exist('default.mat','file')  %If a default pref file exists, load it
        try
            load_pref(handles,1);
        catch                       %If load fails, continue anyway
            %Do nothing on catch (Guide saves a default within the .fig file)
        end
    end

    set(0,'Units','pixels'); %Ensure 'pixel' is the default unit
    guidata(hObject, handles); %save the handles

    clc;        %Clear the command console

    v = version;  %Check version, warn if untested version
    if str2double(v(1:3)) < 7.3
        warning('GUI:VersionCheck','This program has not been tested in MATLAB versions below 7.3.0');
    end

    clear mex;  %Clear any linked mex files
    clear all;  %Clear any workspace variables (not including appdata)

% --- Outputs from this function are returned to the command line.
function varargout = gui_cp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
    varargout{1} = handles.output;

% --- Executes when exit is pressed on the file menu.
function Exit_Callback(hObject, eventdata, handles)
    close('all','hidden');     %Close everything, including any stray hidden windows
    clear all;                 %Clear all the workspace variables
    clear mex;                 %Clear all mex variables

% --- Executes on button press in im_file_browse.
function im_file_browse_Callback(hObject, eventdata, handles)
    %Note: This method is used instead of "'multiselect','on'" as there is a known
    %issue that causes some versions of matlab to stop acknowledging files beyond ~50
    %when this is active
    cur = pwd; %Find the current directory
    if isappdata(0, 'Last_Directory')
        last_dir = getappdata(0, 'Last_Directory');
        if last_dir ~= 0
            if exist(last_dir,'dir')
                cd(last_dir);
            end
        end
    end

    [im_file im_path] = uigetfile({'*.tif;*.tiff;*.bmp;*.jpg;*.jpeg;*.png;*.gif,Image Files (*.tif, *.bmp, *.jpg, *.png, *.gif)'},'Select a frame:');

    if im_file~=0   %Ensure the path is not empty before setting the edit box
        set(handles.im_file_edit, 'String', [im_path,im_file]);
    end

    if im_path ~=0
        if exist(im_path,'dir')
            setappdata(0, 'Last_Directory', im_path);
        end
    end
    cd(cur); %Reset active directory

% --- Executes on button press in out_dir_browse.
function out_dir_browse_Callback(hObject, eventdata, handles)
    cur = pwd; %Find the current directory
    if isappdata(0, 'Last_Directory')
        last_dir = getappdata(0, 'Last_Directory');
        if last_dir ~= 0
            if exist(last_dir,'dir')
                cd(last_dir);
            end
        end
    end

    out_path = uigetdir();                              %Display a UI browse box

    if out_path~=0                                      %Set only if path is valid
        if out_path(length(out_path))~='\'
            set(handles.out_dir_edit, 'String', [out_path,'\']);
        else
            set(handles.out_dir_edit, 'String', out_path);
        end
    end

    if out_path ~=0
        if exist(out_path,'dir')
            setappdata(0, 'Last_Directory', out_path);
        end
    end
    cd(cur); %Reset active directory

% --- Executes on button press in construct_final (Localize Molecules).
function construct_final_Callback(hObject, eventdata, handles)
    try
        if ~get(handles.rolling_ball_cb,'Value')
            if get(handles.custom_mean_bkg, 'Value');   %If saved mbkg is checked, look for a recent save or ask to be pointed to a proper mat file
                if isappdata(0,'last_mbkg_save')
                    last_mbkg = getappdata(0,'last_mbkg_save');
                    answer = questdlg([sprintf('Use most recent Widefield Sum & Mean Background mat file created this session?\n\n'),sprintf('Most Recent:\n'),last_mbkg],' ');
                        if strcmp(answer,'Yes')
                            filename = getappdata(0,'last_mbkg_save');
                        elseif strcmp(answer,'No')
                            cur = pwd; %Find the current directory
                            if isappdata(0, 'Last_Directory')
                                last_dir = getappdata(0, 'Last_Directory');
                                if last_dir ~= 0
                                    if exist(last_dir,'dir')
                                        cd(last_dir);
                                    end
                                end
                            end

                            [fname,pname] = uigetfile({'*.mat','*.mat Files'},'Select Widefield Sum & Mean Background mat file:');

                            if pname ~=0
                                if exist(pname,'dir')
                                    setappdata(0, 'Last_Directory', pname);
                                end
                            end
                            cd(cur); %Reset active directory

                            if fname~=0                     %load the appropiate mat file if valid
                               filename = [pname,fname];
                            else
                               return
                            end
                        else
                            return
                        end

                    setappdata(0,'last_mbkg_save',filename); %set the filename/path in appdata
                    load(filename,'image_sum','field0','n_field');
                else
                    cur = pwd; %Find the current directory
                    if isappdata(0, 'Last_Directory')
                        last_dir = getappdata(0, 'Last_Directory');
                        if last_dir ~= 0
                            if exist(last_dir,'dir')
                                cd(last_dir);
                            end
                        end
                    end

                    [fname,pname] = uigetfile({'*.mat','*.mat Files'},'Select Widefield Sum & Mean Background mat file:');

                    if pname ~=0
                        if exist(pname,'dir')
                            setappdata(0, 'Last_Directory', pname);
                        end
                    end
                    cd(cur); %Reset active directory

                    if fname~=0                     %load the appropiate mat file if valid
                       filename = [pname,fname];
                    else
                       return
                    end

                    setappdata(0,'last_mbkg_save',filename);    %Record the last file used
                    load(filename,'image_sum','field0','n_field');
                end
            else
                [image_sum, field0, n_field] = sum_wf_meanbkg(handles);
            end
        end

        screen = get(0,'ScreenSize');       %Get screensize
        s_width=screen(3);
        s_height=screen(4);
        preview = get(handles.show_preview,'Value');
                
        if get(handles.custom_roi,'Value')
            x_offset   = str2double(get(handles.x_off_edit,'String'));
            y_offset   = str2double(get(handles.y_off_edit,'String'));
            x_size     = str2double(get(handles.x_size_edit,'String'));
            y_size     = str2double(get(handles.y_size_edit,'String'));
        end

        if preview
            %Set up a figure to be used in einzelreader & wf&merge
            h = figure('Position',[10 .13*s_height-64 .93*s_width-10 .91*s_height-35],'Name', ... 
                       'FPALM Analysis','DeleteFcn',@delete_analysis,'UserData',1,'KeyPressFcn',{@change_opt,guidata(gcbo)});
            if ~get(handles.rolling_ball_cb,'Value')
                position3 = [.01,.70,.25,.25];
                plot3 = axes('Position', position3);
                setappdata(plot3,'pos',position3);
                imagesc(image_sum,'HitTest','off')
                if get(handles.custom_roi,'Value')
                    line([x_offset x_offset+x_size x_offset+x_size x_offset x_offset],[y_offset y_offset y_offset+y_size y_offset+y_size y_offset],'Color','r');
                    text(x_offset,y_offset-2,'ROI','Color','r');
                end
                colormap gray
                axis image;
                title(['Widefield Sum: Frames: ',num2str(field0),'-',num2str(n_field)])
                set(plot3,'ButtonDownFcn',{@focus_swap,guidata(gcbo)})
                drawnow;
            end

            if get(handles.adv_cb,'Value')
                einzelreader(handles,h); %use einzelreader with advanced thresholding
            else
                einzelreader2(handles,h); %use einzelreader with simplified thresholding
            end
        else
            if get(handles.adv_cb,'Value')
                einzelreader(handles,0); %use einzelreader with advanced thresholding
            else
                einzelreader2(handles,0); %use einzelreader with simplified thresholding
            end
        end
    catch   %Clean up variables and appdata on error
        clear mex;

        if isappdata(0, 'temp_runtime')
            trt = getappdata(0, 'temp_runtime');
            if exist(trt,'file')
                delete(trt);
            end
            rmappdata(0,'temp_runtime')
        end

        err = lasterror;
        if isempty(findstr('Invalid handle',err.message)) %Ignore 'invalid handle' errors as these are likely the result of intentional manipulation
            rethrow_gui_err(err);
        else
            return
        end
    end

%Change preview options during execution and while figure is selected
function change_opt(hObject, eventdata, handles)
    if eventdata.Key == 's' & strcmp('shift',eventdata.Modifier)
        x = inputdlg(sprintf('Enter new scaling maximum in pixel value:\n(auto = 0)'),'Change scaling limit:',1,{'0'});
        if isempty(x) | isnan(str2double(x))
            %Do nothing
        elseif str2double(x) == 0
            set(handles.mnl_scale_check,'Value',0);
            mnl_scale_check_Callback(0,0,handles);
        else
            set(handles.mnl_scale_check,'Value',1);
            mnl_scale_check_Callback(0,0,handles);
            set(handles.mnl_scale_edit,'String',x);
        end
    elseif eventdata.Key == 'd' & strcmp('shift',eventdata.Modifier)
        x = inputdlg(sprintf('Enter new frame delay in milliseconds (ms):\n(disable = 0)'),'Change frame delay:',1,{'0'});
        if isempty(x) | isnan(str2double(x))
            %Do nothing
        elseif str2double(x) == 0
            set(handles.frame_delay_check,'Value',0);
            frame_delay_check_Callback(0,0,handles);
        else
            set(handles.frame_delay_check,'Value',1);
            frame_delay_check_Callback(0,0,handles);
            set(handles.frame_delay_edit,'String',x);
        end
    elseif eventdata.Key == 'b' & strcmp('shift',eventdata.Modifier)
        if isappdata(0,'focused')
            focused = getappdata(0,'focused');
            if get(handles.show_colorbar,'Value');
                set(handles.show_colorbar,'Value',0);
                axes(focused);
                colorbar('delete');
            else
                set(handles.show_colorbar,'Value',1);
                pos = get(focused, 'Position');
                axes(focused);
                colorbar('Position',[.947, pos(2), .015, pos(4)],'DrawMode','fast')
            end
        end
    end

% --- Executes on button press in compute_mean_bkg_sum.
function compute_mean_bkg_sum_Callback(hObject, eventdata, handles)
    try
        [image field0 n_field meanbkg_all] = sum_wf_meanbkg(handles);   %Obtain mean background image sum

        figure('Name','Widefield Image Sum');             %Create a figure and plot it
        imagesc(image)
        colormap gray
        axis image       %Eliminate white spaceon axis and set to equal
        title(['Widefield Image Sum: Frames: ',num2str(field0),'-',num2str(n_field)])
        drawnow;

        if ~get(handles.custom_mean_bkg,'Value')
            set(handles.custom_mean_bkg,'Value',1);
        end

    catch   %Clean up variables on error
        clear mex;
        err = lasterror;
        if isempty(findstr('Invalid handle',err.message))      %Ignore 'invalid handle' errors as these are most likely the result of intentional manipulation
            rethrow_gui_err(err);
        else
            return
        end
    end

% --- Executes on button press in plot_hist.
function plot_hist_Callback(hObject, eventdata, handles)
    try
        histogram(handles);       %Plot a histogram using plot_histogram.m
    catch
        rethrow_gui_err();
    end

%Executes when figure h is deleted
function delete_analysis(hObject, eventdata)
    if isappdata(gcf,'PrBar_Handle')
        prbar = getappdata(gcf,'PrBar_Handle');

        if ishandle(prbar)  %Remove the progress bar associated with the figure so that stray progress
            delete(prbar);  %bars are not floating around with no associated function or figure
        end

        clear mex;
    end

function render_fpalm_image_Callback(hObject, eventdata, handles)
    try
        answer = questdlg(sprintf('IMPORTANT NOTE: Some variables will be taken from the GUI edit boxes and not the mat file. Output will be saved to the current output directory specified.\n\nContinue?'),'Reconstruct image from mat file:');
        if ~strcmp(answer,'Yes')
            return
        end

        pause(.1);  %Pause to avoid dialog box flickering

        cur = pwd; %Find the current directory
        if isappdata(0, 'Last_Directory')
            last_dir = getappdata(0, 'Last_Directory');
            if last_dir ~= 0
                if exist(last_dir,'dir')
                    cd(last_dir);
                end
            end
        end

        [fname,pname] = uigetfile({'*.mat','*.mat Files'},'Select mat file:');

        if pname ~=0
            if exist(pname,'dir')
                setappdata(0, 'Last_Directory', pname);
            end
        end
        cd(cur); %Reset active directory

        if ~fname
            return
        end

        einzelmat = [pname,fname];

        %Setup appropriate appdata and run pframeit
        setappdata(0,'last_mbkg_save',einzelmat);
        setappdata(0,'temp_runtime',einzelmat);
        setappdata(0,'from_mat_file',1);

        fpalm_render(handles, 0);

        %Clean up variables
        rmappdata(0,'temp_runtime');
        rmappdata(0,'from_mat_file');
        rmappdata(0,'last_mbkg_save');
    catch
        if isappdata(0,'from_mat_file')   %Clean up app data on error
            rmappdata(0,'from_mat_file');
        end
        if isappdata(0,'Save')
            rmappdata(0,'Save');
        end
        if isappdata(0,'temp_runtime')
            rmappdata(0,'temp_runtime');
        end
        if isappdata(0,'last_mbkg_save');
            rmappdata(0,'last_mbkg_save');
        end

        err = lasterror;
        if isempty(findstr('Invalid handle',err.message)) %Ignore handle errors, they are intentional
            rethrow_gui_err(err);
        else
            return
        end
    end
    
% --- Executes on button press in wf_sum_button.
function wf_sum_button_Callback(hObject, eventdata, handles)
    try
        answer = questdlg(sprintf('IMPORTANT NOTE: Some variables will be taken from the GUI edit boxes and not the mat file. Output will be saved to the current output directory specified.\n\nContinue?'),'Reconstruct image from mat file:');
        if ~strcmp(answer,'Yes')
            return
        end

        pause(.1);  %Pause to avoid dialog box flickering

        cur = pwd; %Find the current directory
        if isappdata(0, 'Last_Directory')
            last_dir = getappdata(0, 'Last_Directory');
            if last_dir ~= 0
                if exist(last_dir,'dir')
                    cd(last_dir);
                end
            end
        end

        [fname,pname] = uigetfile({'*.mat','*.mat Files'},'Select mat file:');

        if pname ~=0
            if exist(pname,'dir')
                setappdata(0, 'Last_Directory', pname);
            end
        end
        cd(cur); %Reset active directory

        if ~fname
            return
        end

        einzelmat = [pname,fname];

        sum_wf(handles,einzelmat);

    catch   %Clean up variables on error
        clear mex;
        err = lasterror;
        if isempty(findstr('Invalid handle',err.message))      %Ignore 'invalid handle' errors as these are most likely the result of intentional manipulation
            rethrow_gui_err(err);
        else
            return
        end
    end

% --- Executes on button press in density_plot_button.
function density_plot_button_Callback(hObject, eventdata, handles)
    try
        answer = questdlg(sprintf('IMPORTANT NOTE: Some variables will be taken from the GUI edit boxes and not the mat file. Output will be saved to the current output directory specified.\n\nContinue?'),'Reconstruct image from mat file:');
        if ~strcmp(answer,'Yes')
            return
        end

        pause(.1);  %Pause to avoid dialog box flickering

        cur = pwd; %Find the current directory
        if isappdata(0, 'Last_Directory')
            last_dir = getappdata(0, 'Last_Directory');
            if last_dir ~= 0
                if exist(last_dir,'dir')
                    cd(last_dir);
                end
            end
        end

        [fname,pname] = uigetfile({'*.mat','*.mat Files'},'Select mat file:');

        if pname ~=0
            if exist(pname,'dir')
                setappdata(0, 'Last_Directory', pname);
            end
        end
        cd(cur); %Reset active directory

        if ~fname
            return
        end

        einzelmat = [pname,fname];

        density_plot(handles,einzelmat);

    catch   %Clean up variables on error
        clear mex;
        err = lasterror;
        if isempty(findstr('Invalid handle',err.message))      %Ignore 'invalid handle' errors as these are most likely the result of intentional manipulation
            rethrow_gui_err(err);
        else
            return
        end
    end

%Executes on Load Preferences menu selection
function load_pref_Callback(hObject, eventdata, handles)
    try
        load_pref(handles,0);
    catch
        rethrow_gui_err();
    end

%Executes on Save Preferences menu selection
function save_pref_Callback(hObject, eventdata, handles)
    try
        save_pref(handles,0);
    catch
        rethrow_gui_err();
    end

% --- Executes on About button press in help menu.
function About_Callback(hObject, eventdata, handles)
    msgbox([sprintf('Report bugs, errors, corrections, & suggestions to:\n\nmathew.parent@umit.maine.edu\ntravis.gould@umit.maine.edu\nsam.t.hess@umit.maine.edu\n\n'),'Last Modified: 07/7/08'],'About','Help')

% ------Executes on Documentation button press in help menu.
function GUI_doc_Callback(hObject, eventdata, handles)
    winopen('GUI_instr070708.pdf');
    
    % --- Executes on button press in set_tol_button.
function set_tol_button_Callback(hObject, eventdata, handles)
    set(handles.render_uipanel,'Visible','off');
    set(handles.tol_uipanel,'Visible','on');
    
% --- Executes on button press in tol_ok.
function tol_ok_Callback(hObject, eventdata, handles)
    set(handles.tol_uipanel,'Visible','off');
    set(handles.render_uipanel,'Visible','on');
    
function save_as_default_Callback(hObject, eventdata, handles)
    save_pref(handles,1);   %Save current state, 1 signals these should be saved as default

%Load preferences from a saved file
%NOTE: This function can also attempt to load preferences from an
%      einzelreader mat save.
function load_pref(handles,default)
    if default
        filename = 'default.mat';
    else
        cur = pwd; %Find the current directory
        if isappdata(0, 'Last_Directory')
            last_dir = getappdata(0, 'Last_Directory');
            if last_dir ~= 0
                if exist(last_dir,'dir')
                    cd(last_dir);
                end
            end
        end

        [fname, pname] = uigetfile({'*.mat','*.mat Files'},'Load Preferences: Enter filename and location:');

        if pname ~=0
            if exist(pname,'dir')
                setappdata(0, 'Last_Directory', pname);
            end
        end
        cd(cur); %Reset active directory

        if ~fname
            return
        else
            filename = [pname,fname];
        end
    end

    load(filename);

    if exist('savegui','var')   %Check for the dummy variable
        set(handles.custom_roi,'Value',custom_roi);
            custom_roi_Callback(0, 0, handles);
        set(handles.show_preview,'Value',show_preview);
            show_preview_Callback(0, 0, handles);    
        set(handles.custom_frame,'Value',custom_frame);
            custom_frame_Callback(0, 0, handles);    
        set(handles.custom_mean_bkg,'Value',custom_mean_bkg);
            custom_mean_bkg_Callback(0, 0, handles);      
        set(handles.adv_cb,'Value',use_adv);
            adv_cb_Callback(0, 0, handles);
        if exist('use_lp','var') %for backward compatibility
            set(handles.lp_N_cb,'Value',use_lp);
                lp_N_cb_Callback(0, 0, handles);
        end
        set(handles.write_image_checkbox,'Value',write_to_file);         
        set(handles.rolling_ball_cb,'Value',use_rb);
            rolling_ball_cb_Callback(0, 0, handles);
        if savegui > 2
            set(handles.pixel_photon, 'Value', pixel_photon);
                pixel_photon_Callback(0, 1, handles);
            set(handles.oc_hist,'Value',oc_hist);
            set(handles.oc_wf,'Value',oc_wf);
        end
    elseif exist('total_molecules','var')   %Check for a variable that would be present in einzelreader mat save
        set(handles.custom_roi,'Value',1);
            custom_roi_Callback(0, 0, handles);
        set(handles.show_preview,'Value',1);
            show_preview_Callback(0, 0, handles);    
        set(handles.custom_frame,'Value',1);
            custom_frame_Callback(0, 0, handles);    
        set(handles.custom_mean_bkg,'Value',0);
            custom_mean_bkg_Callback(0, 0, handles);
        set(handles.pixel_photon, 'Value', 0);
            pixel_photon_Callback(0, 1, handles);
    else
        error('Error: Cannot find preferences in file specified.')
    end

    %These variables not required to load, existence check added for backward
    %compatability in mat files that did not store this data
    varlist = {'cam_pix',                    'zero_level',               'exf',...
               'wvlnth',                     'bkgn',                     'pix_to_pho',...
               'iprod_thresh'                'threshold'                 'upper_threshold',...
               'max_pixels_above_threshold', 'n_bright_pixel_threshold', 'box_overlap_factor',...
               'x_offset',                   'y_offset',                 'x_size',...
               'y_size',                     'n_start',                  'n_end',...
               'outpath',                    'ishift',                   'jshift',...
               'NA',                         'rbox',                     'xw_render',...
               'yw_render',                  'size_fac',                 'maxmol',...
               'um_per_pixel',               'psf_scale',                'rb_radius',...
               'FWHM',                       'r0_tol_min',               'r0_tol_max',...
               'lp_tol_min',                 'lp_tol_max',               'N_tol_min',...
               'N_tol_max',                  'frac_unc',                 'max_unc'};
    hndlist = {handles.cam_pix_edit,         handles.zero_lvl,           handles.exf_edit,...
               handles.wvlnth,               handles.bkgn_noise,         handles.pix_to_pho,...
               handles.iprod_thresh_edit,    handles.threshold_edit,     handles.upperthresh_edit,...
               handles.max_above_edit,       handles.min_above_edit,     handles.box_overlap_edit,...
               handles.x_off_edit,           handles.y_off_edit,         handles.x_size_edit,...
               handles.y_size_edit,          handles.from_edit,          handles.to_edit,...
               handles.out_dir_edit,         handles.ishift,             handles.jshift,...
               handles.NA,                   handles.rbox_edit,          handles.xw_render_edit,...
               handles.yw_render_edit,       handles.size_fac_edit,      handles.maxmol_edit,...
               handles.um_per_pixel_edit,    handles.psf_scale_edit      handles.rb_radius_edit,...
               handles.FWHM_edit,            handles.r0_tol_min,         handles.r0_tol_max,...
               handles.lp_tol_min,           handles.lp_tol_max,         handles.N_tol_min,...
               handles.N_tol_max,            handles.frac_unc,           handles.max_unc};
           
    lnth = length(varlist);

    for i=1:lnth
        if exist(varlist{i},'var')
            set(hndlist{i}, 'String', eval(varlist{i}));
        else
            set(hndlist{i}, 'String', 'n/a');
        end
    end
    psf_scale_edit_Callback(0, 0, handles);

%******Special cases******
    if exist('file_info','var')
        if isfield(file_info,'prec')
            set(handles.prec_edit, 'String', file_info.prec);
        end
    end

    if exist('bkg_pc','var')    %For backward compatability
        set(handles.bkg_percent,'String',bkg_pc);
        bkg_pc = str2double(get(handles.bkg_percent,'String'));
        if bkg_pc < 1
            set(handles.bkg_percent,'String',num2str(bkg_pc*100));
        end
    else
        set(handles.bkg_percent,'String','n/a');
    end

    if exist('fs0','var') %Always have a default value for fs
        set(handles.fs,'String',fs0);
    else
        set(handles.fs,'String',50);
    end

    if exist('prec_edit','var')
        if ~isempty(prec_edit)
            set(handles.prec_edit,'String',prec_edit);
        end
    end

    if exist('imagefile','var') %2 ways to load this var
        set(handles.im_file_edit,'String',imagefile);
    elseif exist('infile','var')
        set(handles.im_file_edit,'String',infile);
    end

    if exist('x0_bar','var') && exist('y0_bar','var')
        x0_bar=sbloc(1); %coordinates of scale bar
        y0_bar=sbloc(2);
        barl = [num2str(x0_bar),',',num2str(y0_bar)];
        set(handles.sb_loc,'String',barl); %str2num must be used here
        if exist('sb','var')    %Special formatting for sbloc
            set(handles.sb,'Value',sb);
            sb_Callback(0, 0, handles);            
        end
    end
    
    clear all;

%Save all current preferences to a file
function save_pref(handles,default)
    if default && exist('default.mat','file')
        answer = questdlg('Default preferences already exist. Overwrite?','');

        if ~strcmp(answer,'Yes')
            return
        else
            filename = 'default.mat';
        end
    elseif default
        filename = 'default.mat';
    else
        cur = pwd; %Find the current directory
        if isappdata(0, 'Last_Directory')
            last_dir = getappdata(0, 'Last_Directory');
            if last_dir ~= 0
                if exist(last_dir,'dir')
                    cd(last_dir);
                end
            end
        end
        
        [fname, pname, filter] = uiputfile('*.mat','Save Preferences: Enter filename and location:');
        
        if pname ~=0
            if exist(pname,'dir')
                setappdata(0, 'Last_Directory', pname);
            end
        end
        cd(cur); %Reset active directory

        if fname==0     %If invalid filename, return
            return
        end
        filename = [pname,fname];
    end

    %Save everything typed into the edit boxes, individual variables
    %instead of a save structure used for compatability
    imagefile  = get(handles.im_file_edit,'String');
    outpath    = get(handles.out_dir_edit,'String');
    iprod_thresh = get(handles.iprod_thresh_edit,'String');
    threshold  = get(handles.threshold_edit,'String');
    upper_threshold            = get(handles.upperthresh_edit,'String');
    max_pixels_above_threshold = get(handles.max_above_edit,'String');
    n_bright_pixel_threshold   = get(handles.min_above_edit,'String');
    box_overlap_factor         = get(handles.box_overlap_edit,'String');
    cam_pix    = get(handles.cam_pix_edit,'String');
    x_offset   = get(handles.x_off_edit,'String');
    y_offset   = get(handles.y_off_edit,'String');
    x_size    = get(handles.x_size_edit,'String');
    y_size    = get(handles.y_size_edit,'String');
    bkg_pc     = get(handles.bkg_percent,'String');
    zero_level  = get(handles.zero_lvl,'String');
    exf        = get(handles.exf_edit,'String');
    size_fac   = get(handles.size_fac_edit,'String');
    fs0        = get(handles.fs,'String');
    n_start    = get(handles.from_edit,'String');
    n_end      = get(handles.to_edit,'String');
    file_info.prec  = get(handles.prec_edit,'String');
    pix_to_pho = get(handles.pix_to_pho,'String');
    rbox = get(handles.rbox_edit,'String');
    maxmol = get(handles.maxmol_edit,'String');
    um_per_pixel = get(handles.um_per_pixel_edit,'String');
    psf_scale = get(handles.psf_scale_edit,'String');
    rb_radius = get(handles.rb_radius_edit,'String');
    FWHM = get(handles.FWHM_edit,'String');
    wvlnth = get(handles.wvlnth,'String');
    NA     = get(handles.NA,'String');
    sbloc  = str2num(get(handles.sb_loc,'String')); %str2num must be used here
    x0_bar = sbloc(1); %coordinates of scale bar
    y0_bar = sbloc(2);
    bkgn   = get(handles.bkgn_noise,'String');
    jshift = get(handles.jshift,'String');
    ishift = get(handles.ishift,'String');
    xw_render = get(handles.xw_render_edit,'String');
    yw_render = get(handles.yw_render_edit,'String');
    r0_tol_min = get(handles.r0_tol_min,'String');
    r0_tol_max = get(handles.r0_tol_max,'String');
    lp_tol_min = get(handles.lp_tol_min,'String');
    lp_tol_max = get(handles.lp_tol_max,'String');
    N_tol_min = get(handles.N_tol_min,'String');
    N_tol_max = get(handles.N_tol_max,'String');
    frac_unc = get(handles.frac_unc,'String');
    max_unc = get(handles.max_unc,'String');

    %Save the state of the checkboxes
    custom_roi   = get(handles.custom_roi,'Value');
    show_preview = get(handles.show_preview,'Value');
    custom_frame = get(handles.custom_frame,'Value');
    custom_mean_bkg = get(handles.custom_mean_bkg,'Value');
    pixel_photon = get(handles.pixel_photon,'Value');
    sb      = get(handles.sb,'Value');
    oc_hist = get(handles.oc_hist,'Value');
    oc_wf   = get(handles.oc_wf,'Value');
    write_to_file = get(handles.write_image_checkbox,'Value');
    use_rb   = get(handles.rolling_ball_cb,'Value');
    use_adv   = get(handles.adv_cb,'Value');
    use_lp   = get(handles.lp_N_cb,'Value');

    %Dummy variable to identify save file type
    savegui = 3;

    %Make sure that handles are not saved by parsing through a list of all
    %variables with regular expressions, and removing those that match
    list = who;
    match = regexp(list,'^plot|handles|\<h\>|\<w\>');
    for i=length(match):-1:1
        if match{i}
            list(i)=[];
        end
    end

    %Save the data to a mat file
    save(filename,list{:})

    if default
       msgbox('Default preferences have been saved ("\default.mat") and will be loaded automatically on program start.','Default preferences saved.') 
    end

    %Clear everything else
    clear all;

% Executes on 'show region' button press
function showregion_Callback(hObject, eventdata, handles)
    try
        file     = get(handles.im_file_edit,'String');
        y_offset = str2double(get(handles.y_off_edit,'String')); %Load info from edit boxes
        x_offset = str2double(get(handles.x_off_edit,'String'));
        x_size   = str2double(get(handles.x_size_edit,'String'));
        y_size   = str2double(get(handles.y_size_edit,'String'));
        ishift   = str2double(get(handles.ishift,'String'));
        jshift   = str2double(get(handles.jshift,'String'));

        image0   = zeros(y_size,x_size);

        if get(handles.custom_frame,'Value') %Preview image using "from:" edit box if possible
                                         %this allows the user to choose which frame is previewed
            try %Attempt to read the file as if it were multi-page
                index = str2double(get(handles.from_edit,'String'));
                if index == 1
                    image1 = imread(file,2);
                end
                image1 = imread(file,index);
            catch
                try %Try reading as a series file if multipage fails
                    file_info = get_file_info2(file, handles, 1);
                    image1 = imread([file_info.path,'\',file_info.part_name,num2str(file_info.start,file_info.prec),file_info.ext]);
                catch %Try reading the file with no index(1st frame) on errors
                    image1 = imread(file);
                end
            end
        else
            image1 = imread(file);
        end

        image0(1:y_size,1:x_size) = image1(y_offset:y_offset+y_size-1,x_offset:x_offset+x_size-1);
        clear image1;
    
        axes(handles.roi_pre);
        imagesc(image0)
        colormap gray
        axis('off','image')
        drawnow
    catch
        rethrow_gui_err();
    end

function combine_Callback(hObject, eventdata, handles)
    v = version;  %Check version, warn if untested version
    if str2double(v(1:3)) < 7.3
        warning('GUI:VersionCheck','This function uses "mutliselect" and may truncate file lists > 50 files or cause errors in MATLAB versions below 7.3.0');
    end

    try
        cur = pwd; %Find the current directory
        if isappdata(0, 'Last_Directory')
            last_dir = getappdata(0, 'Last_Directory');
            if last_dir ~= 0
                if exist(last_dir,'dir')
                    cd(last_dir);
                end
            end
        end

        [im_list im_path] = uigetfile({'*.tif;*.tiff;*.bmp;*.png, Image Files'},'Select frames to combine:','Multiselect','on');

        if im_path ~=0
            if exist(im_path,'dir')
                setappdata(0, 'Last_Directory', im_path);
            end
        end
        cd(cur); %Reset active directory

        if numel(im_list) > 1 && iscellstr(im_list) %If cancel button pressed or only 1 file selected
                                                    %then there is no need to continue
            file_info = get_file_info2([im_path,im_list{1}], handles, 1);

            answer = inputdlg('Enter filename and extension to save the combined series to:','Save to:',1,{[file_info.name,'_catenated',file_info.ext]});
            if isempty(answer)
                return
            end

            out_file = [im_path,char(answer)];

            if exist(out_file,'file')   %If file already exists, prompt for overwrite
                answer = questdlg('File already exists. Overwrite?','');

                if ~strcmp(answer,'Yes')
                    return
                else
                    delete(out_file);
                end
            end

            combine(im_path, out_file, im_list);
        end
    catch
        err = lasterror;
        if isempty(findstr('Invalid handle',err.message)) %Ignore handle errors, they are intentional
            rethrow_gui_err(err);
        else
            return
        end
    end

function auto_bkgn_noise_Callback(hObject, eventdata, handles)
    try
        file = get(handles.im_file_edit,'String');
        ppp = str2double(get(handles.pix_to_pho,'String'));

        file_info = get_file_info2(file, handles, 1);

        x = file_info.width;
        y = file_info.height;

        if strcmp(file_info.type,'multipage')
            infile = file;
            image0 = imread(infile,file_info.start);
        elseif strcmp(file_info.type,'singlepage')
            index = num2str(file_info.start,file_info.prec);
            infile = [file_info.path,'\',file_info.part_name,index,file_info.ext];
            if ~exist(infile,'file')
                infile = file;
            end
            image0 = imread(infile);
        end

        h = figure('WindowStyle','modal','NumberTitle','off','Name',...
            'Select Region to use in Background Noise Calculation:',...
            'Pointer','cross','Units','Normalized','OuterPosition',[0 0 1 1]);
        set(h,'Units','Pixels');
        imagesc(image0,'HitTest','off');
        axis image
        colormap gray
        drawnow

        [path,name] = fileparts(infile);
        if strcmp(file_info.type,'multipage')
            title(gca,['*\',name,' -- Index:',num2str(file_info.start)],'Interpreter','none');
        else
            title(gca,['*\',name],'Interpreter','none');
        end

        try
            waitforbuttonpress
        catch
            return
        end

        p1 = get(gca,'CurrentPoint');
        rbbox;
        p2 = get(gca,'CurrentPoint');

        if ishandle(h)
            delete(h);
        end

        p1 = double(uint16(round(p1(1,1:2)))+1); %Image indexes start at 1 but CurrentPoint starts at 0,
        p2 = double(uint16(round(p2(1,1:2)))+1); %offset by 1 is to resolve this

        if(p1(1) < p2(1))
            x1 = p1(1);
            x2 = p2(1);
        else
            x2 = p1(1);
            x1 = p2(1);
        end
        if(p1(2) < p2(2))
            y1 = p1(2);
            y2 = p2(2);
        else
            y2 = p1(2);
            y1 = p2(2);
        end

        if x1 < 1
            x1 = 1;
        end
        if y1 < 1
            y1 = 1;
        end
        if x2 > x
            x2 = x;
        end
        if y2 > y
            y2 = y;
        end

        if x2-x1 == 0 && y2-y1 ==0
            return
        end

        i = double(image0(y1:y2,x1:x2))/ppp;
        
        set(handles.bkgn_noise,'String', num2str(std(i(:)),'%4.4g'));

    catch
        rethrow_gui_err();
    end

% Auto Background level
function auto_bkglvl_Callback(hObject, eventdata, handles)
    try
        imagefile = get(handles.im_file_edit,'String');
        file_info = get_file_info2(imagefile, handles, 1);

        if strcmp(file_info.type,'multipage')
            maximum = histogram(handles, imagefile, file_info.bitdepth, file_info.start);
        elseif strcmp(file_info.type,'singlepage')
            index = num2str(file_info.start,file_info.prec);
            infile = [file_info.path,'\',file_info.part_name,index,file_info.ext];
            if ~exist(infile,'file')
                infile = imagefile;
            end
            maximum = histogram(handles, infile, file_info.bitdepth);
        end

        set(handles.zero_lvl,'String',num2str(maximum,'%4.4g'));
    catch
        rethrow_gui_err();
    end

function region_select_Callback(hObject, eventdata, handles)
    try
        file = get(handles.im_file_edit,'String');

        file_info = get_file_info2(file, handles, 1);

        x = file_info.width;
        y = file_info.height;

        if strcmp(file_info.type,'multipage')
            image0 = imread(file,file_info.start);
        elseif strcmp(file_info.type,'singlepage')
            index = num2str(file_info.start,file_info.prec);
            infile = [file_info.path,'\',file_info.part_name,index,file_info.ext];
            if ~exist(infile,'file')
                infile = file;
            end
            image0 = imread(infile);
        end

        h = figure('WindowStyle','modal','NumberTitle','off','Name',...
            'Select Region of Interest:',...
            'Pointer','cross','Units','Normalized','OuterPosition',[0 0 1 1]);
        set(h,'Units','Pixels');
        imagesc(image0,'HitTest','off');
        axis image
        colormap gray
        drawnow

        try
            waitforbuttonpress
        catch
            return
        end

        p1 = get(gca,'CurrentPoint');
        rbbox;
        p2 = get(gca,'CurrentPoint');

        if ishandle(h)
            delete(h);
        end

        p1 = double(uint16(round(p1(1,1:2)))+1); %Image indexes start at 1 but CurrentPoint starts at 0,
        p2 = double(uint16(round(p2(1,1:2)))+1); %offset by 1 is to resolve this

        x1 = p1(1);
        x2 = p2(1);
        y1 = p1(2);
        y2 = p2(2);

        if x1 > x
            x1 = x;
        end
        if y1 > y
            y1 = y;
        end
        if x2 > x
            x2 = x;
        end
        if y2 > y
            y2 = y;
        end

        x_size = abs(x1-x2)+1;
        y_size = abs(y1-y2)+1;
        

        if x_size == 0 && y_size == 0
            return
        end

        if (x2 > x1)
            x_off = x1;
        else
            x_off = x2;
        end

        if (y2>y1)
            y_off = y1;
        else
            y_off = y2;
        end

        set(handles.custom_roi,'Value',1);
        custom_roi_Callback(0,0,handles);
        set(handles.x_size_edit,'String',num2str(x_size));
        set(handles.y_size_edit,'String',num2str(y_size));
        set(handles.x_off_edit,'String',x_off);
        set(handles.y_off_edit,'String',y_off);
        showregion_Callback(0,0,handles);
    catch
        rethrow_gui_err();
    end

function NA_Callback(hObject, eventdata, handles)
    NA=str2double(get(handles.NA,'String'));
    wvlnth=str2double(get(handles.wvlnth,'String'))/1000; %Conversion to um
    set(handles.oneovere2,'String',num2str(0.55*wvlnth/NA/1.17,'%4.4g'));
    cam_pix_size=str2double(get(handles.cam_pix_edit,'String'));
    set(handles.oneovere2_pix,'String',num2str(0.55*wvlnth/NA/1.17/cam_pix_size,'%4.4g'));
    psf_scale=str2double(get(handles.psf_scale_edit,'String'));
    set(handles.psf_scaled_um,'String',num2str(psf_scale*0.55*wvlnth/NA/1.17,'%4.4g'));
    set(handles.psf_scaled_pix,'String',num2str(psf_scale*0.55*wvlnth/NA/1.17/cam_pix_size,'%4.4g'));

function wvlnth_Callback(hObject, eventdata, handles)
    NA=str2double(get(handles.NA,'String'));
    wvlnth=str2double(get(handles.wvlnth,'String'))/1000; %Conversion to um
    set(handles.oneovere2,'String',num2str(0.55*wvlnth/NA/1.17,'%4.4g'));
    cam_pix_size=str2double(get(handles.cam_pix_edit,'String'));
    set(handles.oneovere2_pix,'String',num2str(0.55*wvlnth/NA/1.17/cam_pix_size,'%4.4g'));
    psf_scale=str2double(get(handles.psf_scale_edit,'String'));
    set(handles.psf_scaled_um,'String',num2str(psf_scale*0.55*wvlnth/NA/1.17,'%4.4g'));
    set(handles.psf_scaled_pix,'String',num2str(psf_scale*0.55*wvlnth/NA/1.17/cam_pix_size,'%4.4g'));

function cam_pix_edit_Callback(hObject, eventdata, handles)
    NA=str2double(get(handles.NA,'String'));
    wvlnth=str2double(get(handles.wvlnth,'String'))/1000; %Conversion to um
    set(handles.oneovere2,'String',num2str(0.55*wvlnth/NA/1.17,'%4.4g'));
    cam_pix_size=str2double(get(handles.cam_pix_edit,'String'));
    set(handles.oneovere2_pix,'String',num2str(0.55*wvlnth/NA/1.17/cam_pix_size,'%4.4g'));
    psf_scale=str2double(get(handles.psf_scale_edit,'String'));
    set(handles.psf_scaled_um,'String',num2str(psf_scale*0.55*wvlnth/NA/1.17,'%4.4g'));
    set(handles.psf_scaled_pix,'String',num2str(psf_scale*0.55*wvlnth/NA/1.17/cam_pix_size,'%4.4g'));
    
function psf_scale_edit_Callback(hObject, eventdata, handles)
    NA=str2double(get(handles.NA,'String'));
    wvlnth=str2double(get(handles.wvlnth,'String'))/1000; %Conversion to um
    set(handles.oneovere2,'String',num2str(0.55*wvlnth/NA/1.17,'%4.4g'));
    psf_scale=str2double(get(handles.psf_scale_edit,'String'));
    set(handles.psf_scaled_um,'String',num2str(psf_scale*0.55*wvlnth/NA/1.17,'%4.4g'));
    cam_pix_size=str2double(get(handles.cam_pix_edit,'String'));
    set(handles.oneovere2_pix,'String',num2str(0.55*wvlnth/NA/1.17/cam_pix_size,'%4.4g'));
    set(handles.psf_scaled_pix,'String',num2str(psf_scale*0.55*wvlnth/NA/1.17/cam_pix_size,'%4.4g'));
    
function iprod_thresh_edit_Callback(hObject, eventdata, handles)
    iprod = str2double(get(hObject, 'String'));
    set(handles.threshold_edit,'String',round(iprod*.7));
    set(handles.upperthresh_edit,'String',iprod);

function sb_Callback(hObject, eventdata, handles)
    if get(handles.sb,'Value')
        set(handles.sb_loc,'Enable','on');
    else
        set(handles.sb_loc,'Enable','off');
    end
    
function mnl_scale_check_Callback(hObject, eventdata, handles)
    if get(handles.mnl_scale_check,'Value')
        set(handles.mnl_scale_edit,'Enable','on');
    else
        set(handles.mnl_scale_edit,'Enable','off');
    end

function frame_delay_check_Callback(hObject, eventdata, handles)
    if get(handles.frame_delay_check,'Value')
        set(handles.frame_delay_edit,'Enable','on');
    else
        set(handles.frame_delay_edit,'Enable','off');
    end

function show_colorbar_Callback(hObject, eventdata, handles)
    if isappdata(0,'focused')
        focused = getappdata(0,'focused');
        if ishandle(focused)
            if ~get(handles.show_colorbar,'Value');
                set(hObject,'Value',0);
                axes(focused);
                colorbar('delete');
            else
                set(hObject,'Value',1);
                pos = get(focused, 'Position');
                axes(focused);
                colorbar('Position',[.947, pos(2), .015, pos(4)],'DrawMode','fast')
            end
        end
    end

function localization_toggle_Callback(hObject, eventdata, handles)
    set(handles.localization_toggle,'color',[0.753, 0.753 0.753]);
    set(handles.background_toggle,'color',[0.941,0.941,0.941]);
    set(handles.roi_toggle,'color',[0.941,0.941,0.941]);
    set(handles.optics_toggle,'color',[0.941,0.941,0.941]);
    set(handles.preview_toggle_tab,'color',[0.941,0.941,0.941]);
    set(handles.render_toggle,'color',[0.941,0.941,0.941]);

    set(handles.localization_uipanel,'Visible','on');
    set(handles.rio_uipanel,'Visible','off');
    set(handles.background_uipanel,'Visible','off');
    set(handles.optics_uipanel,'Visible','off');
    set(handles.preview_uipanel,'Visible','off');
    set(handles.render_uipanel,'Visible','off');
    set(handles.tol_uipanel,'Visible','off');
    
function background_toggle_Callback(hObject, eventdata, handles)
    set(handles.localization_toggle,'color',[0.941,0.941,0.941]);
    set(handles.background_toggle,'color',[0.753, 0.753 0.753]);
    set(handles.roi_toggle,'color',[0.941,0.941,0.941]);
    set(handles.optics_toggle,'color',[0.941,0.941,0.941]);
    set(handles.preview_toggle_tab,'color',[0.941,0.941,0.941]);
    set(handles.render_toggle,'color',[0.941,0.941,0.941]);

    set(handles.localization_uipanel,'Visible','off');
    set(handles.rio_uipanel,'Visible','off');
    set(handles.background_uipanel,'Visible','on');
    set(handles.optics_uipanel,'Visible','off');
    set(handles.preview_uipanel,'Visible','off');
    set(handles.render_uipanel,'Visible','off');
    set(handles.tol_uipanel,'Visible','off');
    
function roi_toggle_Callback(hObject, eventdata, handles)
    set(handles.localization_toggle,'color',[0.941,0.941,0.941]);
    set(handles.background_toggle,'color',[0.941,0.941,0.941]);
    set(handles.roi_toggle,'color',[0.753, 0.753 0.753]);
    set(handles.optics_toggle,'color',[0.941,0.941,0.941]);
    set(handles.preview_toggle_tab,'color',[0.941,0.941,0.941]);
    set(handles.render_toggle,'color',[0.941,0.941,0.941]);

    set(handles.localization_uipanel,'Visible','off');
    set(handles.rio_uipanel,'Visible','on');
    set(handles.background_uipanel,'Visible','off');
    set(handles.optics_uipanel,'Visible','off');
    set(handles.preview_uipanel,'Visible','off');
    set(handles.render_uipanel,'Visible','off');
    set(handles.tol_uipanel,'Visible','off');
    
function optics_toggle_Callback(hObject, eventdata, handles)
    set(handles.localization_toggle,'color',[0.941,0.941,0.941]);
    set(handles.background_toggle,'color',[0.941,0.941,0.941]);
    set(handles.roi_toggle,'color',[0.941,0.941,0.941]);
    set(handles.optics_toggle,'color',[0.753, 0.753 0.753]);
    set(handles.preview_toggle_tab,'color',[0.941,0.941,0.941]);
    set(handles.render_toggle,'color',[0.941,0.941,0.941]);

    set(handles.localization_uipanel,'Visible','off');
    set(handles.rio_uipanel,'Visible','off');
    set(handles.background_uipanel,'Visible','off');
    set(handles.optics_uipanel,'Visible','on');
    set(handles.preview_uipanel,'Visible','off');
    set(handles.render_uipanel,'Visible','off');
    set(handles.tol_uipanel,'Visible','off');

function preview_toggle_tab_Callback(hObject, eventdata, handles)
    set(handles.localization_toggle,'color',[0.941,0.941,0.941]);
    set(handles.background_toggle,'color',[0.941,0.941,0.941]);
    set(handles.roi_toggle,'color',[0.941,0.941,0.941]);
    set(handles.optics_toggle,'color',[0.941,0.941,0.941]);
    set(handles.preview_toggle_tab,'color',[0.753, 0.753 0.753]);
    set(handles.render_toggle,'color',[0.941,0.941,0.941]);
   
    set(handles.localization_uipanel,'Visible','off');
    set(handles.rio_uipanel,'Visible','off');
    set(handles.background_uipanel,'Visible','off');
    set(handles.optics_uipanel,'Visible','off');
    set(handles.preview_uipanel,'Visible','on');
    set(handles.render_uipanel,'Visible','off');
    set(handles.tol_uipanel,'Visible','off');
    
function render_toggle_Callback(hObject, eventdata, handles)
    set(handles.localization_toggle,'color',[0.941,0.941,0.941]);
    set(handles.background_toggle,'color',[0.941,0.941,0.941]);
    set(handles.roi_toggle,'color',[0.941,0.941,0.941]);
    set(handles.optics_toggle,'color',[0.941,0.941,0.941]);
    set(handles.preview_toggle_tab,'color',[0.941,0.941,0.941]);
    set(handles.render_toggle,'color',[0.753, 0.753 0.753]);
   
    set(handles.localization_uipanel,'Visible','off');
    set(handles.rio_uipanel,'Visible','off');
    set(handles.background_uipanel,'Visible','off');
    set(handles.optics_uipanel,'Visible','off');
    set(handles.preview_uipanel,'Visible','off');
    set(handles.render_uipanel,'Visible','on');
    set(handles.tol_uipanel,'Visible','off');
    
function pixel_photon_Callback(hObject, eventdata, handles)
    if get(handles.pixel_photon,'Value')
%         set(handles.pix_to_pho,'Enable','on');
        set(handles.pix1,'String','photons','HitTest','off');
        set(handles.pix2,'String','photons','HitTest','off');
        set(handles.pix3,'String','photons','HitTest','off');

        iprod = str2double(get(handles.iprod_thresh_edit,'String'));
        thresh = str2double(get(handles.threshold_edit,'String'));
        upperthresh = str2double(get(handles.upperthresh_edit,'String'));
        conv_factor = str2double(get(handles.pix_to_pho,'String'));

        if ~isnan(conv_factor) && isempty(eventdata)
            set(handles.iprod_thresh_edit,'String',num2str(round(iprod/conv_factor)));
            set(handles.threshold_edit,'String',num2str(round(thresh/conv_factor)));
            set(handles.upperthresh_edit,'String',num2str(round(upperthresh/conv_factor)));
        end
    else
%         set(handles.pix_to_pho,'Enable','off');
        set(handles.pix1,'String','pixel value','HitTest','off');
        set(handles.pix2,'String','pixel value','HitTest','off');
        set(handles.pix3,'String','pixel value','HitTest','off');

        iprod = str2double(get(handles.iprod_thresh_edit,'String'));
        thresh = str2double(get(handles.threshold_edit,'String'));
        upperthresh = str2double(get(handles.upperthresh_edit,'String'));
        conv_factor = str2double(get(handles.pix_to_pho,'String'));

        if ~isnan(conv_factor) && isempty(eventdata)
            set(handles.iprod_thresh_edit,'String',num2str(round(iprod*conv_factor)));
            set(handles.threshold_edit,'String',num2str(round(thresh*conv_factor)));
            set(handles.upperthresh_edit,'String',num2str(round(upperthresh*conv_factor)));
        end    
    end

% --- Executes on button press in custom_roi.
function custom_roi_Callback(hObject, eventdata, handles)
    if get(handles.custom_roi,'Value')  %Enable/Disable appropriate edit boxes on ROI checkbox change
        set(handles.y_off_edit, 'Enable', 'on');
        set(handles.x_off_edit, 'Enable', 'on');
        set(handles.x_size_edit, 'Enable', 'on');
        set(handles.y_size_edit, 'Enable', 'on');
        set(handles.showregion,'Enable','on');
    else
        set(handles.y_off_edit, 'Enable', 'off');
        set(handles.x_off_edit, 'Enable', 'off');
        set(handles.x_size_edit, 'Enable', 'off');
        set(handles.y_size_edit, 'Enable', 'off');        
        set(handles.showregion,'Enable','off');
    end

% --- Executes on button press in show_preview.
function show_preview_Callback(hObject, eventdata, handles)
    if get(handles.show_preview,'Value')     %Enable/Disable edit box associated with preview checkbox
        set(handles.fs,'Enable','on');
    else
        set(handles.fs,'Enable','off');
    end

% --- Executes on button press in custom_frame.
function custom_frame_Callback(hObject, eventdata, handles)
    if get(handles.custom_frame,'Value')
        set(handles.from_edit,'Enable','on');
        set(handles.to_edit,'Enable','on');
    else
        set(handles.from_edit,'Enable','off');
        set(handles.to_edit,'Enable','off');
    end
    
    % --- Executes on button press in rolling_ball_cb.
function rolling_ball_cb_Callback(hObject, eventdata, handles)
    if get(handles.rolling_ball_cb,'Value')
        set(handles.FWHM_edit,'Enable','on');
        set(handles.rb_radius_edit,'Enable','on');
        set(handles.bkg_percent,'Enable','off');
        set(handles.zero_lvl,'Enable','off');
        set(handles.compute_mean_bkg_sum,'Enable','off');
        set(handles.custom_mean_bkg,'Enable','off');
        set(handles.oc_wf,'Enable','off');
        set(handles.oc_wf,'Value',0);
        set(handles.custom_mean_bkg,'Value',0);
    else
        set(handles.rb_radius_edit,'Enable','off');
        set(handles.FWHM_edit,'Enable','off');
        set(handles.bkg_percent,'Enable','on');
        set(handles.zero_lvl,'Enable','on');
        set(handles.compute_mean_bkg_sum,'Enable','on');
        set(handles.custom_mean_bkg,'Enable','on');
        set(handles.oc_wf,'Enable','on');
    end
    
    % --- Executes on button press in adv_cb (Advanced Thresholding).
function adv_cb_Callback(hObject, eventdata, handles)
    if get(handles.adv_cb,'Value')
        set(handles.threshold_edit,'Enable','on');
        set(handles.upperthresh_edit,'Enable','on');
        set(handles.max_above_edit,'Enable','on');
        set(handles.min_above_edit,'Enable','on');
        set(handles.text14,'Enable','on');
        set(handles.text15,'Enable','on');
        set(handles.text19,'Enable','on');
        set(handles.text20,'Enable','on');
        set(handles.text167,'Enable','on');
        set(handles.text168,'Enable','on');
        set(handles.pix2,'Enable','on');
        set(handles.pix3,'Enable','on');
        set(handles.r0_tol_max,'Enable','off');
        set(handles.r0_tol_min,'Enable','off');
    else
        set(handles.threshold_edit,'Enable','off');
        set(handles.upperthresh_edit,'Enable','off');
        set(handles.max_above_edit,'Enable','off');
        set(handles.min_above_edit,'Enable','off');
        set(handles.text14,'Enable','off');
        set(handles.text15,'Enable','off');
        set(handles.text19,'Enable','off');
        set(handles.text20,'Enable','off');
        set(handles.text167,'Enable','off');
        set(handles.text168,'Enable','off');
        set(handles.pix2,'Enable','off');
        set(handles.pix3,'Enable','off');
        set(handles.r0_tol_max,'Enable','on');
        set(handles.r0_tol_min,'Enable','on');
    end
    
    
% --- Executes on button press in lp_N_cb.
function lp_N_cb_Callback(hObject, eventdata, handles)
    if get(handles.lp_N_cb,'Value')
        set(handles.lp_tol_max,'Enable','on');
        set(handles.lp_tol_min,'Enable','on');
        set(handles.text184,'Enable','on');
        set(handles.text179,'Enable','on');
        set(handles.N_tol_max,'Enable','off');
        set(handles.N_tol_min,'Enable','off');
        set(handles.text183,'Enable','off');
        set(handles.text185,'Enable','off');
    else
        set(handles.lp_tol_max,'Enable','off');
        set(handles.lp_tol_min,'Enable','off');
        set(handles.text184,'Enable','off');
        set(handles.text179,'Enable','off');
        set(handles.N_tol_max,'Enable','on');
        set(handles.N_tol_min,'Enable','on');
        set(handles.text183,'Enable','on');
        set(handles.text185,'Enable','on');
    end

%-------------------------------------------------------------
%   UNUSED CREATEFCNs - These execute when the gui initializes.
%       Most are here to set the proper colors in edit boxes and
%       were automatically created by GUIDE, they are not vital
%       to the operation of this program.
%-------------------------------------------------------------
function fs_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function out_dir_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function zero_lvl_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function bkg_percent_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function box_overlap_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function threshold_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function upperthresh_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function max_above_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function min_above_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function cam_pix_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function iprod_thresh_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function from_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function to_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function pix_x_size_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function exf_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function rk_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function kw_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function x_off_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function y_off_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function x_size_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function im_file_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function prec_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function tab_listbox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function pix_to_pho_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function sb_loc_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
      set(hObject,'BackgroundColor','white');
    end

function bkgn_noise_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function wvlnth_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function jshift_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function ishift_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function NA_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function weight_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function oneovere2_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function mnl_scale_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function frame_delay_edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
       set(hObject,'BackgroundColor','white');
    end

%-------------------------------------------------------------    
%   UNUSED CALLBACKs - These are required by different versions
%       of matlab which somtimes have trouble ignoring missing
%       callbacks.
%-------------------------------------------------------------
function custom_mean_bkg_Callback(hObject, eventdata, handles)

function im_file_edit_Callback(hObject, eventdata, handles)

function x_off_edit_Callback(hObject, eventdata, handles)

function y_off_edit_Callback(hObject, eventdata, handles)

function x_size_edit_Callback(hObject, eventdata, handles)

function box_overlap_edit_Callback(hObject, eventdata, handles)

function min_above_edit_Callback(hObject, eventdata, handles)

function max_above_edit_Callback(hObject, eventdata, handles)

function upperthresh_edit_Callback(hObject, eventdata, handles)

function threshold_edit_Callback(hObject, eventdata, handles)

function bkg_percent_Callback(hObject, eventdata, handles)

function zero_lvl_Callback(hObject, eventdata, handles)

function exf_edit_Callback(hObject, eventdata, handles)

function rk_edit_Callback(hObject, eventdata, handles)

function kw_edit_Callback(hObject, eventdata, handles)

function fs_Callback(hObject, eventdata, handles)

function from_edit_Callback(hObject, eventdata, handles)

function to_edit_Callback(hObject, eventdata, handles)

function out_dir_edit_Callback(hObject, eventdata, handles)

function file_Callback(hObject, eventdata, handles)

function Help_Callback(hObject, eventdata, handles)

function other_options_Callback(hObject, eventdata, handles)

function pix_to_pho_Callback(hObject, eventdata, handles)

function bkgn_noise_Callback(hObject, eventdata, handles)

function sb_loc_Callback(hObject, eventdata, handles)

function oc_hist_Callback(hObject, eventdata, handles)
    
function oc_wf_Callback(hObject, eventdata, handles)

function jshift_Callback(hObject, eventdata, handles)

function ishift_Callback(hObject, eventdata, handles)

function weight_Callback(hObject, eventdata, handles)

function oneovere2_Callback(hObject, eventdata, handles)

function prec_edit_Callback(hObject, eventdata, handles)

function frame_delay_edit_Callback(hObject, eventdata, handles)

function mnl_scale_edit_Callback(hObject, eventdata, handles)

function y_size_edit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function y_size_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit109_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit109_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit110_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit110_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in sb.
function checkbox26_Callback(hObject, eventdata, handles)

function edit111_Callback(hObject, eventdata, handles)

function edit111_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rbox_edit_Callback(hObject, eventdata, handles)

function rbox_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function yw_render_edit_Callback(hObject, eventdata, handles)

function yw_render_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xw_render_edit_Callback(hObject, eventdata, handles)

function xw_render_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function size_fac_edit_Callback(hObject, eventdata, handles)

function size_fac_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxmol_edit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function maxmol_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function um_per_pixel_edit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function um_per_pixel_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function write_image_checkbox_Callback(hObject, eventdata, handles)

function oneovere2_pix_Callback(hObject, eventdata, handles)

function oneovere2_pix_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function psf_scale_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function psf_scaled_um_Callback(hObject, eventdata, handles)

function psf_scaled_um_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function psf_scaled_pix_Callback(hObject, eventdata, handles)

function psf_scaled_pix_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rb_radius_edit_Callback(hObject, eventdata, handles)

function rb_radius_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FWHM_edit_Callback(hObject, eventdata, handles)

function FWHM_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r0_tol_min_Callback(hObject, eventdata, handles)

function r0_tol_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r0_tol_max_Callback(hObject, eventdata, handles)

function r0_tol_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lp_tol_min_Callback(hObject, eventdata, handles)

function lp_tol_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lp_tol_max_Callback(hObject, eventdata, handles)

function lp_tol_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function N_tol_min_Callback(hObject, eventdata, handles)

function N_tol_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function N_tol_max_Callback(hObject, eventdata, handles)

function N_tol_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text187_CreateFcn(hObject, eventdata, handles)

function frac_unc_Callback(hObject, eventdata, handles)

function frac_unc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_unc_Callback(hObject, eventdata, handles)

function max_unc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


