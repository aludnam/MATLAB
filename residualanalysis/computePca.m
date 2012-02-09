function [pcacoef] =computePca(namefile,ncoef)
% [pcacoef] =computePca(namefile,ncoef)

r=load(namefile);
if isfield(r,'peval')
    if isfield(r.peval, 'data_path')
        if isfield(r.peval,'data_file')
            d=load(['~/' r.peval.data_path '/' r.peval.data_file]);
        else 
            d=load(r.peval.data_path); % there was some confusion what is path and waht is file
        end
    else
        error('Can not find the location of the data.')
    end
else 
    d=load('dpixc'); % Data stored directlu in the directory
end

pcacoef = pca(reshape(d.dpixc,r.peval.nx*r.peval.ny,r.peval.nt),ncoef);