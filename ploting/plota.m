function plota(a,peval,savethis)
% plota(a,peval,savethis)

% figure
% plot(a','.:')
% xlabel('time');
% ylabel('a');
% grid on
% setfontsizefigure(12)
% if savethis 
%     SaveImageFULL([peval.res_dir '/a'])
% end

figure
bar(mean(a,2))
hold on
errorbar(1:size(a,1), mean(a,2),std(a,[],2),'+r')
xlabel('component')
ylabel('mean(intensity)')
grid on
setfontsizefigure(12)
if savethis 
    SaveImageFULL([peval.res_dir '/a_mean'])
end
