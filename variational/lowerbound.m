function l=lowerbound(dvec, w, alpha, beta, a, b)
% l=lowerbound(dvec, w, alpha, beta, a, b) (1 x I)
% Computes lower bound for variational approximation (Buntine & Jakulin
% DCA 2006)
% dvec: data - each column is one image (JxI)
% w: basis - each column is one basis vector (component, psf); (JxK)
% alpha, beta: parameters of the Gamma distribution for latent h (1xK each)
% a, b: parameters of the approximative Gamma distribution (KxI each)
% J=#pixels
% K=#components
% I=#images.

% this is the same as in update_n ...
% expectation of <log(l_k)>:
e_l = psi(a) - log(b); % (K x I) 

% n-update:
[ntmp, z]=update_ntmp(w, a, b); % z (Jx1xI) normalization constant from n update

t1=sum(e_l.*bsxfun(@minus, alpha', a),1);   %(1xI)
t2=sum(dvec.*log(squeeze(z)),1);            %(1xI)

% approximation of the log-gamma function
t31=gammalogapprox(a)-a.*log(b);            % (KxI)
t32=gammalogapprox(alpha)-alpha.*log(beta); % (1xK)
t3=sum(bsxfun(@minus, t31, t32'),1);        % (1xI)
% This t4 is not really necessary as it is only data dependent:
t4=sum(factorialapprox(dvec),1);            % (1xI) 
l=t1+t2+t3-t4; %(1xI)