function imstiled (imagein,handlefig,cmap,titleshow,sizevec,steponemore,colorframevec,clim)
% imstiled (imagein,handlefig,cmap,titleshow,sizevec,steponemore,colorframevec))
% Tiles array of images
% handlefig: figure handle, if 0 then handlefig = gcf, if empyt ([]) then
% creates a new figure;
% cmap : colormap
% titleshow : shows title, if set to [] no title is shown
% sizevec : defines size of the subplot. If set to [] then automatic.
% steponemore : (0 default) steps into next window after plotting the last image. 
% colorframevec: a vector indicating which tiles shoudl be framed in
% colored frame: 0-no color, 1-red color, 2-green color, 3-blue color.
% clim: [low high] is the intensity limits of the image. clim=[] fulls the colormap (default)
% (default is no colored frames)
% 
% Example:  im=rand(50,50,6);
%           figure; 
%           imstiled(im,[],'gray',[1:6],[2,3],[],[1 2 3 0 0 0])
% This will create 2x3 tiles from each frame of im, in a gray scale with
% numbers of frames (1:6) in each tile, and making the first frame red,
% second green, third blue adn 4:6 without colored frame. 

if ~exist('steponemore','var');
    steponemore=0;
end

if ~exist('colorframevec','var');
    colorframevec=zeros(size(imagein,3));
end

if ~exist('clim','var')
    clim = []; % image will be scaled to full colormap;
end

if nargin<2 handlefig=0; end
if nargin<3 cmap=[]; end
if nargin<4 titleshow = []; end

d3=size(imagein,3);
d3tmp=d3;
if steponemore
    d3tmp=d3+1;
end

if ~isempty(titleshow)
    if length(titleshow) == 1       
        titlename = num2cell(1:d3); % default title - number of figure
    else %specified ttile as a vector...        
        titlename = num2cell(titleshow);
        titleshow = 1;
    end
end


if nargin<5 % automatic
    sizevec = [];
end
d3=size(imagein,3);
d3tmp=d3;
if steponemore
    d3tmp=d3+1;
end

if isempty(sizevec)
    a = round(sqrt(d3tmp));
    b = ceil(d3tmp/a);
else % defined
    a = sizevec(1);
    b = sizevec(2);
end


if handlefig; 
     figure(handlefig);
else
%      figure(gcf)
end

for ii=1:d3
    subplot(a,b,ii)
    ims(imagein(:,:,ii),cmap,0,1,clim);
    if titleshow
%         title(num2str(titlename{ii}));
%         xlabel(num2str(titlename{ii}),'fontsize',5);
% uncomment this to have it in left top corner:
%         t1 = min(5,floor(size(imagein,1)/5)); %uncomment this to have it in left top corner
% uncomment this to have it in right top corner:
        t1 = min(size(imagein,1)-15,floor(size(imagein,1)-size(imagein,1)/5));
        t2 = min(10,floor(size(imagein,2)/5));
        text(t1, t2, num2str(titlename{ii}),'color','w' ,'fontsize',10)
    end
    if colorframevec(ii)==1
        set(gca,'xcolor','r','linewidth',2)
        set(gca,'ycolor','r','linewidth',2)        
    elseif colorframevec(ii)==2
        set(gca,'xcolor','g','linewidth',2)
        set(gca,'ycolor','g','linewidth',2)        
    elseif colorframevec(ii)==3
        set(gca,'xcolor','b','linewidth',2)
        set(gca,'ycolor','b','linewidth',2)        
    elseif colorframevec(ii)==4
        set(gca,'xcolor','r','linewidth',2)
        set(gca,'ycolor','g','linewidth',2)        
    end 
end

% enable to plot something in another window...
if steponemore 
    subplot(a,b,ii+1)
end