lambda = 655; NA=1.2; pixelsize=106;
raylim=.61*lambda/NA/pixelsize;
% Plotting the variance: 
figure;
clear('l')
hold on
li=length(int1_multi);
intmat=cell2mat(int1_multi);
p=2;
q=2;
plot(repmat(l2(q:end)',1,li-p+1),bsxfun(@times,intmat(p:end)/intmat(p),vardintout(q:end,p:end)));
plot(l2(q:end),2*vard(q:end,p),'--k');
setforsave(gcf,2)
jj=1;
for ii=p:li
    l{jj}=['\alpha=' num2str(1/(int1_multi{ii}/int1_multi{p})) ' (\Lambda=' num2str(int1_multi{ii}) 'photons)'];
    jj=jj+1;
end
l{length(l)+1}='static';
legend(l)
vline2(raylim,'k--')
hold off
xlabel('d [pix]')
ylabel('\alpha \times var(d) [pix^2]')
setfontsizefigure(12)

if savethis 
    SaveImageFULL('images/ComarisonStaticVsBlinking_var')
end

% Plotting the sqrt(variance)

figure; 
clear('l')
hold on
li=length(int1_multi);
intmat=cell2mat(int1_multi);
p=2;
q=5;
plot(repmat(l2(q:end)',1,li-p+1),sqrt(bsxfun(@times,intmat(p:end)/intmat(p),vardintout(q:end,p:end))));
plot(l2(q:end),sqrt(2*vard(q:end,p)),'--k');
setforsave(gcf,2)
jj=1;
for ii=p:li
    l{jj}=['\alpha=' num2str(1/(int1_multi{ii}/int1_multi{p})) ' (\Lambda=' num2str(int1_multi{ii}) 'photons)'];
    jj=jj+1;
end
l{length(l)+1}='static';
legend(l)
vline2(raylim,'k--')
hold off
xlabel('d [pixels]')
ylabel('\alpha \times \surd (var(d)) [pix]')
setfontsizefigure(12)
if savethis 
    SaveImageFULL('images/ComarisonStaticVsBlinking_sqrtvar')
end

% Plotting the components of the Fisher information matrix

figure; 
clear('l')
ristat=reshape(I3d,4,size(I3d,3),4);
ri= reshape(Iintout,4,size(Iintout,3),4);
clear('ccal')
hold on
li=length(int1_multi);
intmat=cell2mat(int1_multi);
p=2;
q=1;
plot(repmat(l2(q:end)',1,li-p+1),bsxfun(@times,1./(intmat(p:end)/intmat(p)),squeeze(ri(1,q:end,p:end))));
plot(l2(q:end),.5*squeeze(ristat(1,q:end,p)),'k');
plot(repmat(l2(q:end)',1,li-p+1),bsxfun(@times,1./(intmat(p:end)/intmat(p)),squeeze(ri(3,q:end,p:end))),'--');
plot(l2(q:end),.5*squeeze(ristat(3,q:end,p)),'--k');
setforsave(gcf,2)
jj=1;
for ii=p:li
    l{jj}=['\alpha=' num2str(int1_multi{ii}/int1_multi{p}) ' (\Lambda=' num2str(int1_multi{ii}) 'photons)'];
    jj=jj+1;
end
l{length(l)+1}='static';
legend(l,'location','best')
vline2(raylim,'k--')
hold off
xlabel('d [pix]')
ylabel('\alpha \times Fisher Information [pix^{-2}]')
setfontsizefigure(12)
if savethis 
    SaveImageFULL('images/ComarisonStaticVsBlinking_FisherInfo')
end
