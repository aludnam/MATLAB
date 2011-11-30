function [pcacoef] =computePca(namefile,ncoef)
% [pcacoef] =computePca(namefile,ncoef)

r=load(namefile);
if isfield(r.res, 'data_file')
    d=load(['~/' r.peval.data_path '/' r.peval.data_file]);
else 
    d=load('dpixc');
end

pcacoef = pca(reshape(d.dpixc,r.peval.nx*r.peval.ny,r.peval.nt),ncoef);