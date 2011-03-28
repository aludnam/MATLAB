function bm=computeblinkmat(blinkmat, Nvec)
% bm=computeblinkmat(blinkmat, Nvec)
% Computes sums of blinkmat for the sources located in each pixel (number
% of sources simulated in each pixel is in vector Nvec)

nc=length(Nvec);
bm = zeros(nc, size(blinkmat,2));

indexstart=1;
for ii=1:nc
    indexstop=indexstart + Nvec(ii) - 1;
    bm(ii,:)=sum(blinkmat(indexstart:indexstop,:),1);
    indexstart=indexstop + 1;
end