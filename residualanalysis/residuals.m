ncvec = [3 4 5 6 8 10];
md = mean(dpixc,3);
for ii= 1:length(ncvec)
    load (['DR2_N5_offset100_iter1_avg1_nc' num2str(ncvec(ii)) '/DR2_N5_offset100_iter1_avg1_nc' num2str(ncvec(ii)) '.mat']);
    % load(['R2/R2_nc' num2str(ncvec(ii))]);
    residual = (dpixc-(reshape(res.w*res.h, peval.nx, peval.ny, peval.nt)));
    mr = mean(abs(residual),3);
    mr2(:,:,ii)=mr;
    [ihs, x_mu, y_mu, sig] = plotreswh(res, peval, dpixc, p, 0, 1, 0, 1, 1);
    ca;
    dipshow(mr2(:,:,ii))
    hold on
    scatter(p.x_vec-0.5, p.y_vec-0.5,200,'xr')
    scatter(x_mu-1, y_mu-1,'b')
    SaveImageFULL(['locmeanAbsResid_nc' num2str(ncvec(ii))])
    
    %     dipshow(mr)
    %     SaveImageFULL(['meanAbsResid_nc' num2str(ncvec(ii))])
    cc = corrcoef(mr,md);
    xc(ii) = cc(2,1);
end