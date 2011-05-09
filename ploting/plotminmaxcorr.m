hf=figure; 
% errorbar(ncvec-1, mean(cpmax),std(cpmax),'o-');
plot(ncvec-1, min(cpmax), 'o-');
% vline2(10,'r-')
xlabel('Number of sources');
ylabel('Maximum correlation in residuals')
setforsave(hf);
if savethis 
    SaveImageFULL('./images/MaxCorrInResid_min')
end
