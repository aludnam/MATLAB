function gd = gradddiv(v, w, h, c, s, nx, ny )
% d = gradddiv(v, w, h, c, s )
% c = ncomp x 2 matrix of centers positions
% s = sigma of gaussian

ncomp = size(c,1);
sv = size(v);

[X, Y] = meshgrid(1:nx, 1:ny);

for ii=1:ncomp
    Xc = X-c(ii,1);
    Yc = Y-c(ii,2);
    t = reshape(-1/s(ii)^2*Xc.*Yc.*exp(-((Xc.^2)+(Yc.^2))/(2*s(ii)^2)), nx*ny,1);
    gdt = (1-v./(w*h)').*repmat(w(:,ii)',sv(1),1).*repmat(t,1, sv(2));
    gd(ii) = sum(gdt(:));
end