function ncomp = estimate_ncpomp(dvec)
fprintf('Repair : ncomp = estimate_ncpomp(dvec)\n')
pca(dvec);

ncomp = 40;
fprintf('Manual: ncomp = %g\n', ncomp)

