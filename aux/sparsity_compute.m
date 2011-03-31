function s = sparsity_compute(M)
% s = sparsity_compute(M)

n=numel(M);
l1 = sum(abs(M(:)));
l2 = sum(M(:).^2);

s = (sqrt(n)-l1/sqrt(l2))/(sqrt(n)-1);