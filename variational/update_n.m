function n=update_n(w, a, b)
% n=update_n(w, a, b)
% Computes update for n_jki (JxKxI) in variational approximation (Buntine & Jakulin
% DCA 2006)
% w: basis - each column is one basis vector (component, psf); (J x K)
% a, b: parameters of the approximative Gamma distribution (KxI each) where
% J=#pixels
% K=#components
% I=#images.
  
% non-normalized n and normalization constant
[ntmp, z]=update_ntmp(w, a, b); % ntmp (JxKxI) , z (Jx1xI)

%normalized such that sum(n,2)=1 for each k
n=bsxfun(@rdivide, ntmp, z); % (J x K x I) 


