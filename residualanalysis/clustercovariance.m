ccd = (corrcoef(dveccr'));  
ccds=squareform(1-ccd);
dendrogram(linkage(ccds), 169);
