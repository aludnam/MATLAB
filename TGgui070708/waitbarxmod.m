function fout = waitbarxmod(x,whichbar, varargin)
%   waitbarxmod.m
%   ----------------------------
%   MODIFIED BY: Matthew Parent
%   LAST UPDATE: 6/26/07
%   REASON:      Coding inefficiencies in original waitbar.m were causing
%                unacceptable slowdowns while executed in 'for' loops.
%   CHANGES:     Waitbar no longer updates on every pass. Waitbar will
%                check for changes of .1% or more before updating. XData
%                for comparison is stored in the appdata of the figure
%                with handle 'f' for fast access.
%
%                **Note that appdata is used in place of the innate
%                UserData field because it is approximately 3-5% faster
%
%                This function uses approximately 12-13% of the CPU time
%                that the built-in waitbar.m uses

%   Original script (waitbar.m) by The Mathworks, Inc.
%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1.23.4.7 $  $Date: 2006/06/27 23:13:17 $

if nargin>=2
    if ischar(whichbar) || iscellstr(whichbar)
        type=2; %we are initializing
        name=whichbar;
    elseif isnumeric(whichbar)
        type=1; %we are updating, given a handle
        f=whichbar;

        if ~ishandle(f) %error out on handle error before attempting to update
            error('Error: Invalid handle.');
        end
    else
        error('MATLAB:waitbar:InvalidInputs', ['Input arguments of type ' class(whichbar) ' not valid.'])
    end
elseif nargin==1
    f = findobj(allchild(0),'flat','Tag','TMWWaitbar');

    if isempty(f)
        type=2;
        name='Waitbar';
    else
        type=1;
        f=f(1);
    end
else
    error('MATLAB:waitbar:InvalidArguments', 'Input arguments not valid.');
end

x = max(0,min(100*x,100));

switch type
    case 1,  % waitbar(x)    update
%         children=allchild(0);
%         if((numel(children)>1) && (children(1)~=f))
%             uistack(f,'top');%Set waitbar on top (tends to cause annoying flickering issues)
%         end

        if nargin>2,
            hAxes = findobj(f,'type','axes'); % Update waitbar title
            hTitle = get(hAxes,'title');
            set(hTitle,'string',varargin{1});
        end

        xc = getappdata(f,'UserData');
        if floor(10*x)==floor(10*xc)%Don't waste CPU updating graphics if updated
            return                  %value will change bar by less than .1%
        else
            setappdata(f,'UserData',x);
        end

        p = findobj(f,'Type','patch');
        l = findobj(f,'Type','line');
        if isempty(f) || isempty(p) || isempty(l),
            error('MATLAB:waitbar:WaitbarHandlesNotFound', 'Couldn''t find waitbar handles.');
        end
        xpatch = get(p,'XData');
        xpatch = [0 x x 0];
        set(p,'XData',xpatch)
        xline = get(l,'XData');
        set(l,'XData',xline);

    case 2,  % waitbar(x,name)  initialize
        vertMargin = 0;
        if nargin > 2,
            % we have optional arguments: property-value pairs
            if rem (nargin, 2 ) ~= 0
                error('MATLAB:waitbar:InvalidOptionalArgsPass',  'Optional initialization arguments must be passed in pairs');
            end
        end

        oldRootUnits = get(0,'Units');

        set(0, 'Units', 'points');
        screenSize = get(0,'ScreenSize');

        axFontSize=get(0,'FactoryAxesFontSize');

        pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');

        width = 360 * pointsPerPixel;
        height = 75 * pointsPerPixel;
        pos = [screenSize(3)/2-width/2 screenSize(4)/2-height/2 width height];

        f = figure(...
            'Units', 'points', ...
            'BusyAction', 'queue', ...
            'Position', pos, ...
            'Resize','off', ...
            'CreateFcn','', ...
            'NumberTitle','off', ...
            'IntegerHandle','off', ...
            'MenuBar', 'none', ...
            'Tag','TMWWaitbar',...
            'Interruptible', 'off', ...
            'DockControls', 'off', ...
            'Visible','off');

        setappdata(f,'UserData',x);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set figure properties as passed to the fcn
        % pay special attention to the 'cancel' request
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        visValue = 'on';
        if nargin > 2,
            propList = varargin(1:2:end);
            valueList = varargin(2:2:end);
            cancelBtnCreated = 0;

            visibleExist = strmatch('vis',lower(propList));
            if ~isempty(visibleExist)
                visValue = valueList{visibleExist};
            end

            for ii = 1:length( propList )
                try
                    if strcmpi(propList{ii}, 'createcancelbtn' ) && ~cancelBtnCreated
                        cancelBtnHeight = 23 * pointsPerPixel;
                        cancelBtnWidth = 60 * pointsPerPixel;
                        newPos = pos;
                        vertMargin = vertMargin + cancelBtnHeight;
                        newPos(4) = newPos(4)+vertMargin;
                        callbackFcn = [valueList{ii}];
                        set( f, 'Position', newPos, 'CloseRequestFcn', callbackFcn );
                        cancelButt = uicontrol('Parent',f, ...
                            'Units','points', ...
                            'Callback',callbackFcn, ...
                            'ButtonDownFcn', callbackFcn, ...
                            'Enable','on', ...
                            'Interruptible','off', ...
                            'Position', [pos(3)-cancelBtnWidth*1.4, 7,  ...
                            cancelBtnWidth, cancelBtnHeight], ...
                            'String','Cancel', ...
                            'Tag','TMWWaitbarCancelButton');
                        cancelBtnCreated = 1;
                    else
                        % simply set the prop/value pair of the figure
                        set( f, propList{ii}, valueList{ii});
                    end
                catch
                    disp ( ['Warning: could not set property ''' propList{ii} ''' with value ''' num2str(valueList{ii}) '''' ] );
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        colormap([]);

        axNorm=[.05 .3 .9 .2];
        axPos=axNorm.*[pos(3:4),pos(3:4)] + [0 vertMargin 0 0];

        h = axes('XLim',[0 100],...
            'YLim',[0 1],...
            'Box','on', ...
            'Units','Points',...
            'FontSize', axFontSize,...
            'Position',axPos,...
            'XTickMode','manual',...
            'YTickMode','manual',...
            'XTick',[],...
            'YTick',[],...
            'XTickLabelMode','manual',...
            'XTickLabel',[],...
            'YTickLabelMode','manual',...
            'YTickLabel',[]);

        tHandle=title(name);
        tHandle=get(h,'title');
        oldTitleUnits=get(tHandle,'Units');
        set(tHandle,...
            'Units',      'points',...
            'String',     name);

        tExtent=get(tHandle,'Extent');
        set(tHandle,'Units',oldTitleUnits);

        titleHeight=tExtent(4)+axPos(2)+axPos(4)+5;
        if titleHeight>pos(4)
            pos(4)=titleHeight;
            pos(2)=screenSize(4)/2-pos(4)/2;
            figPosDirty=true;
        else
            figPosDirty=false;
        end

        if tExtent(3)>pos(3)*1.10;
            pos(3)=min(tExtent(3)*1.10,screenSize(3));
            pos(1)=screenSize(3)/2-pos(3)/2;

            axPos([1,3])=axNorm([1,3])*pos(3);
            set(h,'Position',axPos);

            figPosDirty=true;
        end

        if figPosDirty
            set(f,'Position',pos);
        end

        xpatch = [0 x x 0];
        ypatch = [0 0 1 1];
        xline = [100 0 0 100 100];
        yline = [0 0 1 1 0];

        p = patch(xpatch,ypatch,'b','EdgeColor','b','EraseMode','none');
        l = line(xline,yline,'EraseMode','none');
        set(l,'Color','k');%get(gca,'XColor'));
        set(f,'HandleVisibility','callback','visible', visValue);
        set(0, 'Units', oldRootUnits);
end  % case
drawnow;

if nargout==1,
    fout = f;
end