function plotlowerbound(lb,peval,savethis, xlow)
% plotlowerbound(lb,peval,savethis, xlow)
% xlow : (optional) lower bound on x axes (not to show first few pointss in
% lb...)
if ~exist('xlow', 'var') xlow = 1; end

figure
plot(lb(xlow:end))
xlabel('iteration');
ylabel('lower bound');
grid on
setfontsizefigure(12)
if savethis 
    SaveImageFULL([peval.res_dir '/lowerbound'])
end