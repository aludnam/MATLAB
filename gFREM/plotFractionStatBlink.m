
% fraction= sqrt(vardintout./(2*vard));
fraction= vardintout./(2*vard);
q=3;
f=figure; 
plot(l2(q:end), fraction(q:end,:));

xlabel('d [pix]')
ylabel('var(d)/var^{static}(d)')
for ii=1:length(int2_multi)
    l{ii}=['\Lambda=' num2str(int1_multi{ii}) 'photons'];    
end
legend(l,'location', 'best');
setfontsizefigure(12)
hline(1,'k')
setforsave(f,2);
if savethis 
    SaveImageFULL('images/FractionStatBlink')
end

