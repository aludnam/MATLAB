function plotscatter(h,htrue,savethis,filename)
% plotscatter(h,htrue,savethis)
if ~exist('savethis','var')
    savethis = 0;
end
if ~exist('filename','var')
    filename = 'hvshtrue';
end


minv=min(min(h(:)), min(htrue(:)));
maxv=min(max(h(:)), max(htrue(:)));

figure
hold on
for ii=1:min(size(h))
    scatter(h(ii,:)',htrue(ii,:)','.')
end
plot([minv:maxv],[minv:maxv],'--k','linewidth',2)
xlabel('h')
ylabel('htrue')
set(gca,'dataaspectratio',[1 1 1])
grid on
setfontsizefigure(12)
if savethis
    SaveImageFULL(filename)
end