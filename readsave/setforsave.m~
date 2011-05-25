function setforsave(h)
% setforsave(h)
% Sets paramaters of the graph for saving. (linewidt, fontsize, minorticks)
if ~exist('h', 'var')
    h = gcf;
end
gh=get(h);
ga = get(gh.Children);
figure (h)
setfontsizefigure(12)
set(ga.Children, 'linewidth',1)
set(gh.Children,'XMinorTick','on', 'YMinorTick','on')
grid on