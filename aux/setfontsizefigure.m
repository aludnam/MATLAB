function setfontsizefigure(fontsize)
% setfontsizefigure(fontsize)
% set font size (axes, labels) of the current figure
set(get(gca,'XLabel'),'FontSize',fontsize)
set(get(gca,'YLabel'),'FontSize',fontsize)
set(get(gca,'ZLabel'),'FontSize',fontsize)
set(gca,'FontSize',fontsize)