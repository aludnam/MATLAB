lambda = 655; NA=1.2; pixelsize=106;
raylim=.61*lambda/NA/pixelsize;
% Plotting the variance: 
figure;
clear('l')
hold on
li=length(int1_multi);
intmat=cell2mat(int1_multi);
m=1;
q=2;
plot(repmat(l2(q:end)',1,li-m+1),bsxfun(@times,intmat(m:end)/intmat(m),vardintout(q:end,m:end)));
plot(l2(q:end),2*vard(q:end,m),'--k');
setforsave(gcf,2)
jj=1;
for ii=m:li
    l{jj}=['\alpha=' num2str(1/(int1_multi{ii}/int1_multi{m})) ' (\Lambda=' num2str(int1_multi{ii}) 'photons)'];
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
    name = 'images/ComarisonStaticVsBlinking_var';
    if p.offset > 0
        name = [name '_bg' num2str(p.offset)];
    end
    SaveImageFULL(name)
end

% Plotting the sqrt(variance)

figure; 
clear('l')
hold on
li=length(int1_multi);
intmat=cell2mat(int1_multi);
plot(repmat(l2(q:end)',1,li-m+1),sqrt(bsxfun(@times,intmat(m:end)/intmat(m),vardintout(q:end,m:end))));
plot(l2(q:end),sqrt(2*vard(q:end,m)),'--k');
setforsave(gcf,2)
jj=1;
for ii=m:li
    l{jj}=['\alpha=' num2str(1/(int1_multi{ii}/int1_multi{m})) ' (\Lambda=' num2str(int1_multi{ii}) 'photons)'];
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
    name = 'images/ComarisonStaticVsBlinking_sqrtvar';
    if p.offset > 0
        name = [name '_bg' num2str(p.offset)];
    end
    SaveImageFULL(name)
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
q=1;
plot(repmat(l2(q:end)',1,li-m+1),bsxfun(@times,1./(intmat(m:end)/intmat(m)),squeeze(ri(1,q:end,m:end))));
plot(l2(q:end),.5*squeeze(ristat(1,q:end,m)),'k');
plot(repmat(l2(q:end)',1,li-m+1),bsxfun(@times,1./(intmat(m:end)/intmat(m)),squeeze(ri(3,q:end,m:end))),'--');
plot(l2(q:end),.5*squeeze(ristat(3,q:end,m)),'--k');
setforsave(gcf,2)
jj=1;
for ii=m:li
    l{jj}=['\alpha=' num2str(int1_multi{ii}/int1_multi{m}) ' (\Lambda=' num2str(int1_multi{ii}) 'photons)'];
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
    name = 'images/ComarisonStaticVsBlinking_FisherInfo';
    if p.offset > 0
        name = [name '_bg' num2str(p.offset)];
    end
    SaveImageFULL(name)
end
