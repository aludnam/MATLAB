function ploterr(sep, offset,errNMF,varNMF,savethis)
% function ploterr(sep, offset,errNMF,varNMF,savethis)

errorbar(repmat(sep', 1,length(offset)), errNMF', varNMF','o-')
makelegend(offset, 'NMF offset', 'separation [pixels]', 'err loc [pixels]')
xlim([sep(1)-0.1, sep(end)+0.1])
if exist('savethis','var')
        SaveImageFULL('NMF_err',savethis)
end