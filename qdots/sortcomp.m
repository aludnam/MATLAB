function sc = sortcomp(compin)
%nejak to moc nefunguje....
compd = dip_image(compin);
comp = double(gaussf(compd,1));

amax = squeeze(max(max(abs(comp))));
amean = squeeze(mean(mean((abs(comp)))));
v = (amax-amean)./(amax+amean);
[a, b] = sort(v, 'descend');
sc = comp(:,:,b);

imax = squeeze(max(max(sc)));
imin = squeeze(min(min(sc)));
msub = find(abs(imin)>imax); %giggest nonzero value negative
sc(:,:,msub) = -sc(:,:,msub); %flip negative to positive






