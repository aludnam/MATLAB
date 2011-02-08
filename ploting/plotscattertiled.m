function plotscattertiled(h,htrue,savethis,filename)
% plotscattertiled(h,htrue,savethis,filename)
if ~exist('savethis','var')
    savethis = 0;
end
if ~exist('filename','var')
    filename = 'hvshtruetiled';
end

minv=min(min(h(:)), min(htrue(:)));
maxv=min(max(h(:)), max(htrue(:)));

figure
m=min(size(h));
n=min(size(htrue));
ll=1;
for jj=1:n
    for ii=1:m
        subplot(n,m,ll)
        ll=ll+1;
        scatter(htrue(jj,:)',h(ii,:)','.')
        hold on
        plot([minv:maxv],[minv:maxv],'--k','linewidth',2)
        
        xlabel(['htrue' num2str(jj)])
        ylabel(['h' num2str(ii)])
        set(gca,'dataaspectratio',[1 1 1])
        grid on
        xlim([0,maxv])
        ylim([0,maxv])
        
    end
end


if savethis
    SaveImageFULL(filename)
end