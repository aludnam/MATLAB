% this plots var(d) for 'improper' treatment of the blinking situation (not
% integrating out the intensities)
% results vard from teh function varianceFREM.m

q=1;plot(l2(q:end),(bsxfun(@times, [1 2],vard(q:end,:))))
xlim([.5,8])
xlabel('d [pixels]')
ylabel('var(d) [pixels^2]')
setforsave(gcf,2)
setfontsizefigure(12)
legend({'{(1,1)}', '{(1,1),(1,0),(0,1),(0,0)}'})
if savethis
    SaveImageFULL('images/VarianceDistanceEqualVsBlinkingNotIntOut')
end
% 

% plot(l2, 2*.5*1000*vard, l2, .5*1000*vardintout)
% xlim([.5,8])
% xlabel('d [pixels]')
% ylabel('var(d) [pixels^2]')
% setforsave(gcf,2)
% setfontsizefigure(12)
% legend({'{(1,1)}', '{(1,1),(1,0),(0,1),(0,0)}'})
% if savethis
%     SaveImageFULL('images/VarianceDistanceEqualVsBlinkingIntOut')
% end