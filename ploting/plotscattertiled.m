function correlationcoef=plotscattertiled(h,htrue,savethis,filename,showcorrcoef)
% plotscattertiled(h,htrue,savethis,filename)
if ~exist('savethis','var')
    savethis = 0;
end
if ~exist('filename','var')
    filename = 'hvshtruetiled';
end
if ~exist('showcorrcoef','var')
    showcorrcoef =1;
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
        cc=corrcoef(htrue(jj,:)',h(ii,:)');
        correlationcoef(jj,ii)=cc(1,2);
        hold on
        plot([minv:maxv],[minv:maxv],'--k','linewidth',2)
        
        xlabel(['htrue' num2str(jj)])
        ylabel(['h' num2str(ii)])
        set(gca,'dataaspectratio',[1 1 1])
        grid on
        xlim([0,maxv])
        ylim([0,maxv])
        if showcorrcoef
%             text(round(minv+10), round(maxv-10),['Corr. Coef = ' num2str(correlationcoef(jj,ii))]);
            title(['Corr. Coef = ' num2str(correlationcoef(jj,ii))]);
        end
    end
end


if savethis
    SaveImageFULL(filename)   
end

posx=[.5,1.5];
posy=[1.5,0.5]+0.1;
h=hinton(correlationcoef');
for ii=1:size(correlationcoef,1)
    for jj=1:size(correlationcoef,2)
        text(posx(ii),posy(jj),num2str(correlationcoef(ii,jj)),'fontsize',20)
    end
end

        
if savethis
    set(h, 'inverthardcopy', 'off')
    SaveImageFULL([filename, '_hintoncorrocef'])   
end
