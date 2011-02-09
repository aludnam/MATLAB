function [x_mu, y_mu sig] = localizew(w,peval)
% [x_mu, y_mu sig] = localizew(w,peval)
% Localilzes w (#pixels x #components) from NMF model (V=WH)
% peval : parameters (needed peval.nx, peval.ny)

K=size(w,2);
wr=reshape(w,peval.nx, peval.ny,K); %K can be different from peval.ncomp if background is takes as one component...

for ii=1:K %background is not localized...
    [x_mu(ii), y_mu(ii), sig(ii), differ(ii)] = fitgauss2d(wr(:,:,ii));
end