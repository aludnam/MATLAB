function d = distmat (mat)
% d = distmat (mat)
% compute length of vectors in columns of mat
d = sqrt(sum(mat.^2,1));
