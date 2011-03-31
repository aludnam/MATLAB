function wg = makegauss(c, sig, sizevec)
% wg = makegauss(c, sig, sizevec)
% Creates N images of gaussians centered at c (Nx2 matrix). Each image
% contains one gaussian centered at c(ii,:) and with std sig(ii). Each
% image image will be of the size sizevec. 
    
ng = size(c,1); % Number of gaussians
if numel(sig) == 1 % All will have the same sig
    sig = repmat(sig, 1, ng);
end

[X,Y] = meshgrid(0:sizevec(1)-1,0:sizevec(2)-1);
n0=1./(2*pi*sig.^2);

for ii=1:ng
    wg(:,:,ii)=n0(ii)*exp(-((X-c(ii,1)).^2+(Y-c(ii,2)).^2)/(2*sig(ii)^2));
end
wg=normalize(wg);
