function ploterrpersep(sep, offset,errNMF,varNMF,savethis)
% function ploterrpersep(sep, offset,errNMF,varNMF,savethis)
q=repmat(sep', 1,length(offset));
plot(q, errNMF'./q,'o-')
makelegend(offset, 'NMF offset', 'separation [pixels]', 'err loc/separation')
xlim([sep(1)-0.1, sep(end)+0.1])
if exist('savethis','var')
        SaveImageFULL('NMF_err_per_sep',savethis)
end