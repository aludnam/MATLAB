function [cmax,ll,lb,pcacoef, cmax_artificial, cmax_artificial_var] =computeCorrelLlLb(namefile)
% [cmax,ll,lb,pcacoef, cmax_artificial, cmax_artificial_var] =computeCorrelLlLb(namefile)

r=load(namefile);
if isfield(r,'peval')
    if isfield(r.peval, 'data_path')
        if isfield(r.peval,'data_file')
            d=load([r.peval.data_path '/' r.peval.data_file]);
        else 
            d=load(r.peval.data_path); % there was some confusion what is path and waht is file
        end
    else
        error('Can not find the location of the data.')
    end
else 
    d=load('dpixc'); % Data stored directlu in the directory
end

if isfield(r.res,'a')
    r.res.h=r.res.a;
end
model = r.res.w*r.res.h;
resid = (model - reshape(d.dpixc, r.peval.nx*r.peval.ny, r.peval.nt))./sqrt(model);
c = corrcoef(resid');
cmax = max(c(c<1)); % c<1 restrict to off-diagonal elements

% This computer artificial residuals (model-noise(model,poisson)) which can
% be used as a baseline for the residual analysis (1e12 factor is due to 
% matlab implementation of the function imnoise.m):
n_artif = 10; % number of repetition of the artificial residual simulations
cmax_artificial_n = zeros(1,n_artif);
for ii=1:n_artif
    resid_artificial = (model - 1e12*imnoise(1e-12*model))./sqrt(model);
    c_artificial = corrcoef(resid_artificial');
    cmax_artificial_n(ii) = max(c_artificial(c_artificial<1));
end
cmax_artificial = mean(cmax_artificial_n);
cmax_artificial_var = var(cmax_artificial_n);

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