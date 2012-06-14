%This program localizes single molecules and determines their positions
%by fitting a gaussian to their intensity profile.  A position dependent
%background subtraction is used.

function einzelreader(handles, h)
    global xpix ypix wbox psf_w02;

    version = '6/23/08'; %record version (date) of gui being used

    rbox=str2double(get(handles.rbox_edit,'String'));
    wbox=2*rbox+1;
    [xpix,ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
    
    q       = str2double(get(handles.cam_pix_edit,'String')); 
    wvlnth  = str2double(get(handles.wvlnth,'String'))/1000; %convert wavelength from nm to um       
    NA      = str2double(get(handles.NA,'String'));
    psf_scale = str2double(get(handles.psf_scale_edit,'String'));
    
    psf_w0 = psf_scale*0.55*wvlnth/NA/1.17; % 1/e2 radius of PSF in um, use 1/e2 = FWHM/1.17 from Pawley
                                      % with scale factor (20% for "real objective" due to
                                      % measured PSF from Hess and Webb
                                      % 2002)
    psf_std=psf_w0/2; %standard deviation of psf
    psf_w02=(psf_w0/q)*(psf_w0/q); %square of 1/e^2 radius in pixels
    
    yfit_psf=exp(-2*((xpix).*(xpix)+(ypix).*(ypix))/psf_w02);
    npix=sum(sum(yfit_psf));     % area of molecule in square pixels
    
    total_molecules=0;
    xcm_all=zeros(1,1);
    ycm_all=zeros(1,1);
    framenum_all=zeros(1,1);
    n_fail_a0=0;
    n_fail_outbox=0;

    imagefile = get(handles.im_file_edit, 'String');
    outpath = get(handles.out_dir_edit, 'String');
    if outpath(length(outpath))~='/'
        outpath = [outpath,'/'];
    end
    preview = get(handles.show_preview,'Value');

    fs0 = str2double(get(handles.fs,'String')); %Setup for preview frame skipping
    fs1 = fs0;
    
    if get(handles.rolling_ball_cb,'Value') %if using rolling ball subtraction
        %Obtain info on image type, etc
        file_info = get_file_info2(imagefile, handles);
        
        n_end   = file_info.stop;
        n_start = file_info.start;
        
        if get(handles.custom_roi,'Value')
            x_offset   = str2double(get(handles.x_off_edit,'String'));
            y_offset   = str2double(get(handles.y_off_edit,'String'));
            x_size     = str2double(get(handles.x_size_edit,'String'));
            y_size     = str2double(get(handles.y_size_edit,'String'));
        else
            x_offset = 1;
            y_offset = 1;
            x_size   = file_info.width;
            y_size   = file_info.height;
        end
        
        rball=str2double(get(handles.rb_radius_edit,'String')); %radius of rolling ball
        se = strel('ball',rball,rball,0); %structural element, i.e. rolling ball
        
        FWHM=str2double(get(handles.FWHM_edit,'String'));; %FWHM of gaussian smoothing in pixels
        rk=(FWHM)/sqrt(2*log(2)); %1/e^2 smoothing radius in pixels
        kw=20; %kernal width of smoothing function
        [X,Y]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
        kd=sqrt(X.*X+Y.*Y);
        gs=exp(-2*kd.*kd/(rk*rk));
        gs=gs/sum(sum(gs)); %smoothing function normalized to have area = 1
        
    else %else use integrated widefield sum
        %Load wf sum & mean background mat file
        meanbkg_mat_file = getappdata(0,'last_mbkg_save'); %strcat(dir,'analysis\',date,'/',base_name,'sum_wf_meanbkg.mat');
        load(meanbkg_mat_file,'meanbkg_all','image_sum','x_size','y_size','x_offset','y_offset','field0','n_field','file_info');%,'zero_level'); %mat file from sum_wf_meanbkg.m
        
        n_end   = file_info.stop;
        n_start = file_info.start;
        %give warning if starting/ending frames are not compatible with wf sum & mean background mat file 
        if(n_start<field0)
            warndlg(sprintf('Warning: Starting frame is less than corresponding frame used in Compute WF & Mean Background, need to change starting frame or re-run Compute WF & Mean Background'));
            return  
        end
        if(n_end>n_field)
            warndlg(sprintf('Warning: Ending frame is greater than corresponding frame used in Compute WF & Mean Background, need to change ending frame or re-run Compute WF & Mean Background'));
            return  
        end
        
        bkg_pc = str2double(get(handles.bkg_percent,'String'))/100;
        zero_level = str2double(get(handles.zero_lvl,'String'));
        
        %normalize so that ave pixel value = 1
        image_sum_mean=mean(mean(image_sum(y_offset:y_offset+y_size-1,x_offset:x_offset+x_size-1)));
        image_sum_norm = (image_sum(y_offset:y_offset+y_size-1,x_offset:x_offset+x_size-1))/image_sum_mean;
    end
    
    %Load variables from the GUI edit boxes
    %--------------------------------------
    if get(handles.pixel_photon,'Value')
        pix_to_pho = str2double(get(handles.pix_to_pho,'String'));
    else
        pix_to_pho = 1;
    end

    iprod_thresh    = str2double(get(handles.iprod_thresh_edit,'String'))*pix_to_pho; %initial threshold for a "bright object"
    threshold       = str2double(get(handles.threshold_edit,'String'))*pix_to_pho; %this value to be included
    upper_threshold = str2double(get(handles.upperthresh_edit,'String'))*pix_to_pho; %this value to not be excluded

    n_bright_pixel_threshold   = str2double(get(handles.min_above_edit,'String')); %min # of pixels above    
    box_overlap_factor         = str2double(get(handles.box_overlap_edit,'String')); %if center to center closer, don't include either
    max_pixels_above_threshold = str2double(get(handles.max_above_edit,'String')); %max # of pixels above
    
    %Initialize temporary runtime file for passing variables
    %-------------
    for x=1:99
        temp_runtime = [outpath,'temp_runtime_',num2str(x,'%2.2d'),'.mat'];
        if ~exist(temp_runtime,'file')
            break
        end
    end
    setappdata(0,'temp_runtime',temp_runtime);
    %-------------

    %Initialize waitbar and pause component
    w = waitbarxmod(0,'Executing "einzelreader.m" ...','CreateCancelBtn','delete(gcf)');
    set(w,'Name','Progress Bar');
    uicontrol('Style','pushbutton','Parent',w,'String','Pause','Position',[210,10,60,23], ...
              'UserData',1,'Callback',@pause_gui);
    drawnow;   %Draw the extra button immediately    
    pause(.1); %Pause to ensure window completes drawing
    %keepontop('Progress Bar');

    wb_norm = n_end-n_start;
    if wb_norm == 0 %Avoid divide by zero errors
        wb_norm = 1;
    end

    if h~=0                             %If preview window is active, store the handle of the progress bar
        setappdata(h,'PrBar_Handle',w); %in the appdata of the figure
    end

    for fileloop=n_start:n_end
        if strcmp(file_info.type,'singlepage')
            index = num2str(fileloop,file_info.prec);
%             infile = [file_info.path,'/',file_info.part_name,index,file_info.ext];
            infile = [file_info.path,'/',file_info.part_name,index,file_info.ext];
            i1=double(imread(infile));
        elseif strcmp(file_info.type,'multipage')
            i1=double(imread(imagefile,fileloop));
        end

        compl = (fileloop-n_start)/wb_norm;
        drawnow;
        waitbarxmod(compl,w,'Executing "einzelreader.m" ...'); %Rename the waitbar
        
        if get(handles.rolling_ball_cb,'Value') %if using rolling ball subtraction
            i1_gs = uint16(conv2(i1,gs,'same')); %smoothed original
            rb_im = imopen(i1_gs,se);
            bkg   = zeros(y_size,x_size);
            bkg   = double(rb_im(y_offset:y_offset+y_size-1,x_offset:x_offset+x_size-1));  
        else
            bkg    = zeros(y_size,x_size);
            bkg    = image_sum_norm*bkg_pc*meanbkg_all(fileloop)+zero_level;
        end

        iprod = zeros(y_size,x_size);
        iprod = i1(y_offset:y_offset+y_size-1,x_offset:x_offset+x_size-1)-bkg;
        iprod=iprod.*(iprod>0); % set any negative pixel values to zero

        high_pixel_mask = zeros(y_size,x_size);
        high_pixel_xy   = zeros(10000,3);  % x,y,active
        n_high_pixels   = 0;
        pix_val=iprod_thresh;

        while(pix_val >= iprod_thresh) %continue until reach minimum threshold
            pix_val=0;
            for i=rbox+1:y_size-rbox-1
                for j=rbox+1:x_size-rbox-1
                    if high_pixel_mask(i,j)==0 && iprod(i,j)>pix_val
                        pix_val=iprod(i,j); %find brightest pixel
                        high_pixel_y=i;
                        high_pixel_x=j;
                    end
                end
            end
            if(pix_val < iprod_thresh)
                break
            end
            x0box=high_pixel_x-rbox;
            y0box=high_pixel_y-rbox;
            x1box=high_pixel_x+rbox;
            y1box=high_pixel_y+rbox;

            high_pixel_mask(y0box:y1box,x0box:x1box)=1;
            n_high_pixels=n_high_pixels+1;
            high_pixel_xy(n_high_pixels,1)=high_pixel_x;
            high_pixel_xy(n_high_pixels,2)=high_pixel_y;
            high_pixel_xy(n_high_pixels,3)=1; % active
        end

        drawnow;              %**These added to keep waitbar responsive...
                              % Without these 'drawnow's, the queue gets large and
                              % ignores button presses during execution...
                              %(this function is fast when there is no
                              % actual drawing to do)

        waitbarxmod(compl,w); % <-check waitbar (fast,see waitbarxmod.m)

        for i=1:n_high_pixels
          for j=1:i
            if (i~=j & high_pixel_xy(i,3)==1 & high_pixel_xy(j,3)==1)
              dx=high_pixel_xy(i,1)-high_pixel_xy(j,1);
              dy=high_pixel_xy(i,2)-high_pixel_xy(j,2);
              dmin_nearest_box=sqrt(dx*dx+dy*dy);

              if (dmin_nearest_box<box_overlap_factor*rbox) % these boxes overlap
                if j>i
                  high_pixel_xy(j,3)=-3;       % make one of the boxes inactive
                else
                  high_pixel_xy(i,3)=-3;       % make one of the boxes inactive
                end
              end
            end
          end
        end

        drawnow;
        waitbarxmod(compl,w); %check waitbar

        n_boxes=0;
        boxes_xy=zeros(10000,3);

        for i=1:n_high_pixels
           if high_pixel_xy(i,3)==1
               n_boxes=n_boxes+1;
               boxes_xy(n_boxes,1)=high_pixel_xy(i,1);
               boxes_xy(n_boxes,2)=high_pixel_xy(i,2);
               boxes_xy(n_boxes,3)=1;
           end
        end

        drawnow;
        waitbarxmod(compl,w); %check waitbar

        failed=zeros(n_boxes,1);
        for i=1:n_boxes
           x0_box=boxes_xy(i,1)-rbox;
           y0_box=boxes_xy(i,2)-rbox;
           x1_box=boxes_xy(i,1)+rbox;
           y1_box=boxes_xy(i,2)+rbox;

           grab=iprod(y0_box:y1_box,x0_box:x1_box);
           minbkg=min(min(grab)); 
           grab=grab-minbkg;
           
           xm_sum=0;
           ym_sum=0;
           m_sum=0;
           for x=x0_box:x1_box
             for y=y0_box:y1_box
                xind=floor(x);
                yind=floor(y);
                intens=iprod(yind,xind);
                xm_sum=xm_sum+xind*intens;
                ym_sum=ym_sum+yind*intens;
                m_sum=m_sum+intens;
             end
           end

           x_cm(i)=xm_sum/m_sum;
           y_cm(i)=ym_sum/m_sum;
           
           %re-center box around center of mass
           if round(x_cm(i)) > 1+rbox && round(x_cm(i)) < x_size-rbox
               boxes_xy(i,1)=round(x_cm(i));
           end
           if round(y_cm(i)) > 1+rbox && round(y_cm(i)) < y_size-rbox
               boxes_xy(i,2)=round(y_cm(i));
           end

           x0_box=boxes_xy(i,1)-rbox;          
           y0_box=boxes_xy(i,2)-rbox;        
           x1_box=boxes_xy(i,1)+rbox;
           y1_box=boxes_xy(i,2)+rbox;
                      
           grab=zeros(wbox,wbox);
           grab=iprod(y0_box:y1_box,x0_box:x1_box);
           minbkg=min(min(grab)); 
           grab=grab-minbkg;
           grab_sum(i)=sum(sum(grab));
           image_bright_in_box(:,:,i)=grab;

           xc_box=(x0_box+x1_box)*0.5;
           yc_box=(y0_box+y1_box)*0.5;

           xguess=x_cm(i)-xc_box;
           yguess=y_cm(i)-yc_box;

           for n=1:wbox
                for p=1:wbox
                    k=(n-1)*wbox+p;
                    xymerge(k)=0;
                    zmerge(k)=grab(n,p);
                end
           end

           beta=[xguess,yguess,50];
           [betafit,resid,J,COVB,mse] = nlinfit(xymerge,zmerge,@gaussian_merge,beta);
           ci = nlparci(betafit,resid,'covar',COVB); %calculate error estimates on parameters
           ci_err=(ci(:,2)-ci(:,1))/2;
           xf_err(i)=ci_err(1);
           yf_err(i)=ci_err(2);
           a0_err(i)=ci_err(3);

           zfitmerge=gaussian_merge(betafit,xymerge);

           for n=1:wbox
                for p=1:wbox
                    k=(n-1)*wbox+p;
                    zfit(n,p)=zfitmerge(k);
                end
           end

           yf_box(i)=betafit(2)+yc_box;
           xf_box(i)=betafit(1)+xc_box;
           a0_box(i)=betafit(3);
        
           if(a0_box(i) < 0)
                n_fail_a0=n_fail_a0+1;
                failed(i)=1;
           end
           if xf_box(i) > x1_box || xf_box(i) < x0_box || yf_box(i) > y1_box || yf_box(i) < y0_box
                n_fail_outbox=n_fail_outbox+1;
                failed(i)=1;
           end

           drawnow   %Flush event queue to keep cancel/pause button responsive on every pass.
                     %This is the slowest function by far and drawnow barely changes its exec time

           waitbarxmod(compl,w); %check waitbar
        end

        n_pixels_above_threshold=zeros(n_boxes,1);
        n_pixels_above_upper_threshold=zeros(n_boxes,1);

        for i=1:n_boxes
          if (boxes_xy(i,3)==1)
            n_pixels_above_threshold(i)=0;
            for j=1:2*rbox+1
              for k=1:2*rbox+1
                intens=image_bright_in_box(j,k,i);
                if intens>=threshold  
                  n_pixels_above_threshold(i)=n_pixels_above_threshold(i)+1;
                end
                if intens>=upper_threshold
                  n_pixels_above_upper_threshold(i)=n_pixels_above_upper_threshold(i)+1;
                end
              end
            end
          end
        end

        drawnow;
        waitbarxmod(compl,w); %check waitbar

        n_molecules_found=0;

        for i=1:n_boxes
            if boxes_xy(i,3)==1
                if (n_pixels_above_threshold(i)>=n_bright_pixel_threshold && n_pixels_above_upper_threshold(i)<=max_pixels_above_threshold)
                    n_molecules_found=n_molecules_found+1;
                else
                    if (n_pixels_above_threshold(i)<n_bright_pixel_threshold)
                      boxes_xy(i,3)=-2; %draw red box
                    end
                    if (n_pixels_above_threshold(i)>max_pixels_above_threshold)
                      boxes_xy(i,3)=-1; %draw yellow box
                    end
                end
            end
        end

        drawnow;
        waitbarxmod(compl,w); %check waitbar

        if preview  %If preview active, update images on screen
            set(0,'CurrentFigure',h);
            if ~exist('plot2','var')
                position2 = [.01,.37,.25,.25];
                plot2 = axes('Position', position2);
                setappdata(plot2,'pos',position2);
            else
                set(h,'CurrentAxes',plot2);
            end
            i1display = i1(y_offset:y_offset+y_size-1,x_offset:x_offset+x_size-1);
            imagesc(i1display,'HitTest','off');
            axis image;
            colormap gray;
            title(['Last Frame Processed: ', num2str(fileloop)])
            set(plot2, 'ButtonDownFcn',{@focus_swap,guidata(gcbo)});
            clear i1display;

            if ~exist('plot1','var')
                position1 = [.01,.05,.25,.25];
                plot1 = axes('Position', position1,'Color','k','YTickLabel','','XTickLabel','');
                axis image;
                set(plot1,'Position',position1);
                setappdata(plot1,'pos', position1);
            end
            
            if get(handles.rolling_ball_cb,'Value')
                if ~exist('plot3','var')
                    position3 = [.01,.70,.25,.25];
                    plot3 = axes('Position', position3);
                    setappdata(plot3,'pos',position3);
                else
                    set(h,'CurrentAxes',plot3);
             end
                imagesc(rb_im,'HitTest','off')
                colormap gray, axis image;
                title(['Rolling Ball Profile: ',num2str(fileloop)])
                set(plot3,'ButtonDownFcn',{@focus_swap,guidata(gcbo)});
            end

            if ~exist('plot4','var')
                if x_size < y_size
                    positionf = [.18,.065,.87,.87];
                else
                    positionf = [.32,.3,.6,.65];
                end
                plot4 = axes('Position', positionf);
                setappdata(0,'focused',plot4);
            else
                set(h,'CurrentAxes',plot4);
            end

            if get(handles.mnl_scale_check,'Value')
                clim = str2double(get(handles.mnl_scale_edit,'String'));
                imagesc(iprod,'HitTest','off',[0 clim]);
            else
                imagesc(iprod,'HitTest','off');
            end

            axis image;
            colormap gray;
            title(['Einzel Reader: Frame: ',num2str(fileloop)]);
            set(plot4, 'ButtonDownFcn',{@focus_swap,guidata(gcbo)});

            if (get(handles.show_colorbar,'Value') == 1)
                focused = getappdata(0, 'focused');
                set(h,'CurrentAxes', focused);
                pos = get(focused, 'Position');
                colorbar('Position',[.947, pos(2), .015, pos(4)],'DrawMode','fast')
            end

            set(h,'CurrentAxes', plot4);
        end

        if n_boxes>0

            if preview %Only draw if preview is active
                hold on
                draw_boxes(n_boxes,boxes_xy,rbox);
            end

            for i=1:n_boxes
                if boxes_xy(i,3)==1 && failed(i) ~= 1
                    if preview %Only plot if preview is active
                        %plot(x_cm(i),y_cm(i),'.m','HitTest','off'); %centroid
                        %plot(xf_box(i),yf_box(i),'.b','HitTest','off'); %gaussian fit
                    end
                    total_molecules=total_molecules+1;
                    xcm_all(total_molecules)=x_cm(i);
                    ycm_all(total_molecules)=y_cm(i);
                    xf_all(total_molecules)=xf_box(i);
                    yf_all(total_molecules)=yf_box(i);
                    a0_all(total_molecules)=a0_box(i);
                    grab_sum_all(total_molecules)=grab_sum(i);
                    framenum_all(total_molecules)=fileloop;
                    xf_err_all(total_molecules)=xf_err(i);
                    yf_err_all(total_molecules)=yf_err(i);
                    a0_err_all(total_molecules)=a0_err(i);
                end
            end
            nmol_all(fileloop)=n_boxes;%%

            if preview %Remove the hold
                hold off
            end

            drawnow;              %Draw any boxes and update einzelreader image (flush queue and check waitbar if preview not active)
            waitbarxmod(compl,w); %Update the waitbar
        end

        if get(handles.frame_delay_check,'Value')
            delay = str2double(get(handles.frame_delay_edit,'String'))/1000;
            if delay ~= inf
                pause(delay);
            end
        end

        if preview  %If preview active, compute wf & merge on first, last, and nth frames (always)
            if fileloop == n_end
                setappdata(0,'Save',1); %Trigger image save on last frame
                save(temp_runtime);     %Save variables to a temporary file for function passing
%                 fpalm_render_einzelreader(handles,h);
            elseif fileloop == n_start+fs1-1
                fs1=fs1+fs0;
                save(temp_runtime);
%                 fpalm_render_einzelreader(handles,h);
            elseif fileloop == n_start
                save(temp_runtime);
%                 fpalm_render_einzelreader(handles,h);
            end
        end 
    end 

    if ~preview %If preview is off, continue to wf & merge at end and save the result
        plot1 = 0;
        save(temp_runtime);
        setappdata(0,'Save',1);
        fpalm_render_einzelreader(handles,h);
    else        %Otherwise delete progress bar and handle stored in appdata
        delete(w);
        rmappdata(h,'PrBar_Handle');
    end
    
    if get(handles.oc_wf,'Value') %save widefield image if selected
        if (file_info.part_name(length(file_info.part_name)) == '_') %Avoid saving with repeated seperators                
            wf_file = [outpath,'/',file_info.part_name,'wf.tif'];
        else
            wf_file = [outpath,'/',file_info.part_name,'_wf.tif'];
        end
        imwrite(uint8((image_sum/max(image_sum(:)))*255),wf_file,'Compression','none');
    end

%Save Variables
%----------------
    if (file_info.part_name(length(file_info.part_name)) == '_') %Avoid saving with repeated seperators
        out_file = [outpath,file_info.part_name,num2str(n_start),'-',num2str(n_end),'_t',num2str(iprod_thresh),'-',num2str(threshold),'-',num2str(upper_threshold),'_npt',num2str(n_bright_pixel_threshold),'-',num2str(max_pixels_above_threshold),'.mat'];
    else
        out_file = [outpath,file_info.part_name,'_',num2str(n_start),'-',num2str(n_end),'_t',num2str(iprod_thresh),'-',num2str(threshold),'-',num2str(upper_threshold),'_npt',num2str(n_bright_pixel_threshold),'-',num2str(max_pixels_above_threshold),'.mat'];
    end

    if ~exist(outpath,'dir') %Make sure the directory exists
        answer = questdlg('Output directory could not be accessed. Choose a new save location?','Warning!','Yes','Discard','Yes');

        if strcmp(answer,'Yes')
            outpath = uigetdir();
            if ~outpath
                return
            end
            if outpath(length(outpath))~='/'
                outpath = [outpath,'/'];
            end

            if (file_info.part_name(length(file_info.part_name)) == '_')
                out_file = [outpath,file_info.part_name,num2str(n_start),'-',num2str(n_end),'_t',num2str(iprod_thresh),'-',num2str(threshold),'-',num2str(upper_threshold),'_npt',num2str(n_bright_pixel_threshold),'-',num2str(max_pixels_above_threshold),'.mat'];
            else
                out_file = [outpath,file_info.part_name,'_',num2str(n_start),'-',num2str(n_end),'_t',num2str(iprod_thresh),'-',num2str(threshold),'-',num2str(upper_threshold),'_npt',num2str(n_bright_pixel_threshold),'-',num2str(max_pixels_above_threshold),'.mat'];
            end
        else
            return
        end
    end

    prec_edit     = get(handles.prec_edit,'String');                %<- These to be saved with the
    bkgn          = str2double(get(handles.bkgn_noise,'String'));   %<- mat file for loading purposes    
    ppp           = str2double(get(handles.pix_to_pho,'String'));   %<-
    jshift        = str2double(get(handles.jshift,'String'));       %<-
    ishift        = str2double(get(handles.ishift,'String'));       %<-
    exf           = str2double(get(handles.exf_edit,'String'));     %<-   expansion factor
    
    a0_phot=a0_all/ppp;    % peak amplitude of each molecule in photons
    N=npix*a0_phot;        % number of photons for each molecule

    %localization precision in um
    lp2=((psf_std^2)+(q^2)/12)*1./N+8*pi*(psf_std^4)*(bkgn^2)/(q^2)*1./(N.*N);
    lp=sqrt(lp2);
    
    list = who; %Use regular expressions to ensure handles are not saved
    match = regexp(list,'^plot|handles|\<h\>|\<w\>');
    for i=length(match):-1:1
        if match{i}
            list(i)=[];
        end
    end
    save(out_file,list{:}); %Save variables used in einzelreader.m
    setappdata(0,'last_einzel_save',out_file);  %Save location of einzelreader mat file for other functions to use
%----------------
    %Clean up temp files and appdata
    if exist(temp_runtime,'file')
        delete(temp_runtime)
    end
    if isappdata(0,'temp_runtime')
        rmappdata(0,'temp_runtime')
    end