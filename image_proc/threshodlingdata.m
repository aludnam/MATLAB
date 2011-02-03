maskd= max(threshold(dpixc),[],3);
smd = sum(maskd);
maskdncomp = repmat(maskd,1,1,peval.ncomp);
maskdnt = repmat(maskd,1,1,peval.nt);
d1=(reshape(dpixc(maskdnt),smd,1,peval.nt));
peval.ny=1;
peval.nx=smd;
peval.numpix = peval.nx;

%...computation...

d2=newim(maskdncomp);
d2(maskdncomp) = squeeze(reshape(res.w(:,1:peval.ncomp),smd*peval.ncomp,1));
