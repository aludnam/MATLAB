
% fraction= sqrt(vardintout./(2*vard));
% fraction= vardintout./(2*vard);
fraction = 1./sqrt(vardintout./vard);
q=3;
f=figure;
plot(l2(q:end)*p.pixelsize, fraction(q:end,:));

xlabel('d [nm]')
ylabel('r')
for ii=1:length(int2_multi)
    l{ii}=['\Lambda=' num2str(int1_multi{ii}/2) 'photons'];
end
% legend(l,'location', 'best');
legend(l);
setfontsizefigure(12)
hline(1,'k')
setforsave(f,2);
if savethis
    %     SaveImageFULL('images/FractionStatBlink')
    SaveImageFULL('images/FractionStatBlink_2DPSF_bg100')
    
end

