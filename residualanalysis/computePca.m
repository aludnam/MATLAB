function [pcacoef] =computePca(namefile,ncoef)
% [pcacoef] =computePca(namefile,ncoef)

r=load(namefile);
d=load(['~/' r.peval.data_path '/' r.peval.data_file]);

pcacoef = pca(reshape(d.dpixc,r.peval.nx*r.peval.ny,r.peval.nt),ncoef);