function setforsave(h)
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