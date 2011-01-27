function n=update_n(w, a, b)
% n=update_n(w, a, b)
% Computes update for n_jk in variational approximation (Buntine & Jakulin
% DCA 2006)
% w: basis - each column is one basis vector (component, psf); (#pixels=J x
% #basis=K)
% a, b: parameters of the approximative Gamma distribution (1xK each)

j=size(w,1); % number of pixels
e_l = psi(a) - log(b); % expectation of <log(l_k)> (1xK)
ntmp=w.*repmat(e_l,j,1);
z=sum(ntmp,1); % normalization constant (1 x K)

n=ntmp./repmat(z,j,1); % (J x K)


