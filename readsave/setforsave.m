function setforsave(h,linewidthvalue)
% setforsave(h,linewidthvalue)
% Sets paramaters of the graph for saving. (linewidt, fontsize, minorticks)
if ~exist('linewidthvalue','var')
    linewidthvalue = 1;
end
if ~exist('h', 'var')
    h = gcf;
end
gh=get(h);
for ii=1:length(gh.Children)
    ga{ii} = get(gh.Children(ii));
end
figure (h)
setfontsizefigure(12)
for ii=1:length(gh.Children)
set(ga{ii}.Children, 'linewidth',linewidthvalue)

end
set(gh.Children,'XMinorTick','on', 'YMinorTick','on')
grid on