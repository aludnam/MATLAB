function plota(a,peval,savethis)
% plota(a,peval,savethis)

figure
plot(a','.:')
xlabel('time');
ylabel('a');
grid on
setfontsizefigure(12)
if savethis 
    SaveImageFULL([peval.res_dir '/a'])
end