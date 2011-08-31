function [cmax,ll,lb,pcacoef] =computeCorrelLlLb(namefile)
% [cmax,ll,lb] =computeCorrLlLb(namefile)

r=load(namefile);
d=load(['~/' r.peval.data_path '/' r.peval.data_file]);
if isfield(r.res,'a')
    r.res.h=r.res.a;
end
resid = (r.res.w*r.res.h - reshape(d.dpixc, r.peval.nx*r.peval.ny, r.peval.nt))./sqrt(r.res.w*r.res.h);

c = (corrcoef(resid'));
cmax = max(c(c<1));
% cpmin(jj,ii) = min(cp(cp<1));
% cpmean(jj,ii) = mean(cp(cp<1));
% cpmeanpos(jj,ii) = mean(cp(and(cp>0,cp<1)));
% cpmeanneg(jj,ii) = mean(cp(cp<0));
ll=loglikelihoodPoisson(reshape(d.dpixc,r.peval.nx*r.peval.ny,r.peval.nt),r.res.w*r.res.h);
if isfield(r.res, 'lb')
    lb = r.res.lb;
else
    lb = 0;
end

pcacoef = pca(reshape(d.dpixc,r.peval.nx*r.peval.ny,r.peval.nt),20);