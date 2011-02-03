function plotindividualw(w,peval,savethis)
% plotindividualw(w,peval,savethis)

wr=reshape(w, peval.nx, peval.ny, peval.ncomp);

for ii=1:peval.ncomp
    dipshow(wr(:,:,ii));
    if savethis
        SaveImageFULL([peval.res_dir '/w' num2str(ii)],'e');
    end
end

