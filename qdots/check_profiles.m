% saves images of individual estimated components....
function varout = check_profiles(sep5, offset5, pathname, prename, niter)
p.path = pathname;


for kk5 = 1: length(offset5)
    for ll5=1 : length(sep5)
        p.namedir = [prename num2str(100*sep5(ll5)) 'offset_' num2str(offset5(kk5))];
        cd ([p.path p.namedir])
        load ([p.namedir '-iter_1'])
        %         ims(psf);
        %         SaveImageFULL('psf', 'pf');
        is = imarfun('imsum',dpixc_ind);
        prof = double(is(:,15));
        varout(kk5, ll5,:) = prof;
% % %         for mm5=1:niter
% % %             varout = bg;
% % %             imstiled(icapixNMF{mm5}(:,:,1:2));
% % %             close all
% % %         end
    end
end