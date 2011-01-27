function a=update_a(dvec,n, alpha)
% a=update_a(dvec,n, alpha)
% Computes update for a_k in variational approximation (Buntine & Jakulin
% DCA 2006)
% dvec: data - each column is one image (#pixels=J x #documents=I)
% #basis=K)
% alpha, beta: parameters of the Gamma distribution for latent h (1xK each)

a = sum(sum(n'*dvec,2)' + alpha);