function [x_mu, y_mu sig] = localizewpix(w)
% [x_mu, y_mu sig] = localizewpix(w)
K=size(w,3)
for ii=1:K %background is not localized...
    [x_mu(ii), y_mu(ii), sig(ii), differ(ii)] = fitgauss2d(w(:,:,ii));
end