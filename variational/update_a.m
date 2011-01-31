function a=update_a(dvec,n,alpha)
% a=update_a(dvec,n, alpha)
% Computes update for a_ki (KxI) in variational approximation (Buntine & Jakulin
% DCA 2006)
% dvec: data - each column is one image (J x I)
% n: (J x K x I)
% alpha, beta: parameters of the Gamma distribution for latent h (1xK each)
% J=#pixels
% K=#components
% I=#images.

[j,k,i]=size(n);

dvec_reshaped = reshape(dvec,j,1,i); % (Jx1xI)
atmp = sum(bsxfun(@times,dvec_reshaped,n),1); % (1xKxI)
a=bsxfun(@plus, squeeze(atmp), alpha'); % (KxI)