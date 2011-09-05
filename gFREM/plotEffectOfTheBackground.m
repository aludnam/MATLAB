f=figure; 
semilogx([0:10:100, 200:100:1000],sqrt([vard_bg04, vard_bg])*p.pixelsize,'--k')
hold on
semilogx([0:10:100, 200:100:1000],sqrt([vardintout_bg04, vardintout_bg])*p.pixelsize,'-k')
xlabel('background [photons]')
ylabel('FREM sqrt(var(d)) [nm]')
l{1}='static sources';
l{2}='blinking sources';
legend (l)
setfontsizefigure(12)
setforsave(f,2);
if savethis
    name = ['EffectOfTheBackground_int' num2str(int1_multi{1}/2) '_d' num2str(l2*10)];
    SaveImageFULL(['images/' name])
end