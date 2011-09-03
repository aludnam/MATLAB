% lambda = 655; NA=1.2; pixelsize=106;
raylim=.61*p.lambda/p.NA;
% Plotting the variance: 
bgvec = [0 10 100 1000];
vard = [vard_bg0, vard_bg10, vard_bg100, vard_bg1000];
vardintout = [vardintout_bg0, vardintout_bg10, vardintout_bg100, vardintout_bg1000];
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
for ii=m:length(bgvec)
    % there is a factor of 1/2 for intensity (int1_multi) as this will
    % reflect the mean number of photons....
    l{jj}=['Background' nim2str(bgvec(ii)) 'photons (\Lambda=' num2str(int1_multi/2) 'photons)'];
    jj=jj+1;
end
legend(l)
% vline2(raylim,'k--')
hold off
xlabel('Distance d [nm]')
ylabel('FREM sqrt(var(d)) [nm]')
setfontsizefigure(12)

if savethis 
    name = 'images/ComarisonStaticVsBlinking_PSF2D_FREM_differentBg';    
    SaveImageFULL(name)
end
