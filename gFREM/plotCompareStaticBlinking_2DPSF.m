% lambda = 655; NA=1.2; pixelsize=106;
raylim=.61*p.lambda/p.NA;
% Plotting the variance: 

clear('l')
li=length(int1_multi);
m=1;
q=2;
l2rep_nm = repmat(l2(q:end)',1,li-m+1)*p.pixelsize;
figure;
semilogy(l2rep_nm,p.pixelsize*sqrt(vardintout(q:end,m:end)),'-')
hold on
semilogy(l2rep_nm,p.pixelsize*sqrt(vard(q:end,m:end)),'--')

% compensation for photons numbers included in computation
setforsave(gcf,2)
jj=1;
for ii=m:li
    % there is a factor of 1/2 for intensity (int1_multi) as this will
    % reflect the mean number of photons....
    l{jj}=['\Lambda=' num2str(int1_multi{ii}/2) 'photons'];
    jj=jj+1;
end
legend(l)
% vline2(raylim,'k--')
hold off
xlabel('Distance d [nm]')
ylabel('FREM sqrt(var(d)) [nm]')
setfontsizefigure(12)

if savethis 
    name = 'images/ComarisonStaticVsBlinking_PSF2D_FREM';
    if p.offset > 0
        name = [name '_bg' num2str(p.offset)];
    end
    SaveImageFULL(name)
end

% % Plotting the components of the Fisher information matrix
% 
% figure; 
% clear('l')
% ristat=reshape(I3d,4,size(I3d,3),4);
% ri= reshape(Iintout,4,size(Iintout,3),4);
% clear('ccal')
% li=length(int1_multi);
% intmat=cell2mat(int1_multi);
% q=1;
% l2rep_nm = repmat(l2(q:end)',1,li-m+1)*p.pixelsize;
% semilogy(l2rep_nm,squeeze(ri(1,q:end,m:end)),'-o');
% hold on
% semilogy(l2rep_nm,squeeze(ri(3,q:end,m:end)),'-s');
% semilogy(l2(q:end)*p.pixelsize,squeeze(ristat(1,q:end,m:end)),'--o');
% semilogy(l2(q:end)*p.pixelsize,squeeze(ristat(3,q:end,m:end)),'--s');
% setforsave(gcf,2)
% jj=1;
% for ii=m:li
%     l{jj}=['\Lambda=' num2str(int1_multi{ii}) 'photons'];
%     jj=jj+1;
% end
% l{length(l)+1}='static';
% legend(l,'location','best')
% % vline2(raylim,'k--')
% hold off
% xlabel('d [pix]')
% ylabel('\alpha \times Fisher Information [pix^{-2}]')
% setfontsizefigure(12)
% if savethis 
%     name = 'images/ComarisonStaticVsBlinking_FisherInfo_PSF2D';
%     if p.offset > 0
%         name = [name '_bg' num2str(p.offset)];
%     end
%     SaveImageFULL(name)
% end
