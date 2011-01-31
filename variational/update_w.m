function w=update_w(dvec, n)
% w=update_w(dvevc, n)
% Computes update for w_jk (JxK) in variational approximation (Buntine & Jakulin
% DCA 2006)
% dvec: (J x I) data - each column is one image 
% n: (J x K x I )
% J=#pixels
% K=#components
% I=#images.

[j,k,i]=size(n);

dvec_reshaped = reshape(dvec,j,1,i);
wtmp = squeeze(sum(bsxfun(@times,dvec_reshaped,n),3));
w=bsxfun(@rdivide, wtmp, sum(wtmp,1));







