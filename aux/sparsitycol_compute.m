function s = sparsitycol_compute(M)
% sparsitycol_compute is the Hoyer's sparsity of the columns of M. 
% s = sparsitycol_compute(M)


n=size(M,1);
l1 = sum(abs(M),1);
l2 = sum(M.^2,1);

s = (sqrt(n)-l1./sqrt(l2))/(sqrt(n)-1);