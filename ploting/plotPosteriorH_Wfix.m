% x1=reshape(blinkmat, 1,peval.ncomp*peval.nt);
% x1=Hres;
x2=[p.x_vec, p.y_vec];
nmax=1000;
step=nmax/20;
h1=(0:step:nmax);
h2=(0:step:nmax);
h1(1)=eps;
h2(1)=eps;
p1=peval;
p1.nt=1;

Winit='trueW';
%  Winit='wrongW';
for kk=1:4
    
    d1=dpixc(:,:,kk);
    
    if strcmp(Winit,'trueW')
        x2in = x2;
    elseif strcmp(Winit,'wrongW')
        x2in = Wres;
    else
        error('wrong name')
    end
    for ii=1:length(h1)
        for jj=1:length(h2)
            x1in = log([h1(ii), h2(jj)]);
            P(ii,jj)=loglikGaPExpPriorH ([x1in x2in] ,reshape(d1, p1.nx*p1.ny, p1.nt), p1.sigmaPSF, 1, peval.beta, p1);
        end
    end
    %%
    figure
    [a(kk),b(kk)]=find(P==min(P(:)));
    [C,h] = contour(h1, h2, -P,50);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
    set(gca,'xtick',(0:step:nmax))
    set(gca,'ytick',(0:step:nmax))
    grid on
    hold on
    scatter (h1(b(kk)),h2(a(kk)),'xr')
    scatter (blinkmat(2,kk), blinkmat(1,kk),'ob')
    scatter (blinkmat(1,kk), blinkmat(2,kk),'ob')
    t1=['true h1: ';'true h2: '];
    title([t1, num2str(blinkmat(:,kk))])
    lg{1}='P(D|W,H)';
    lg{2}='max';
    lg{3}='true';
    legend(lg)
    xlabel('h1');
    ylabel('h2');
    
    if savethis
        SaveImageFULL(['posteriorHfixW-Winit' Winit '-timeSlice' num2str(kk)])
    end
end
