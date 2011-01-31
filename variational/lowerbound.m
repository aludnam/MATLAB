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
%
% There should be another term log(dvec!) but it is constant during
% evaluation...

[j,i]=size(dvec);
k=size(w,2);

% this is the same as in update_n ...
% expectation of <log(l_k)>:
e_l = psi(a) - log(b); % (K x I) 
% Expand into another dimension for multiplication with w...
e_lreshaped=reshape(e_l,1,k,i); %  (1 x K x I) 

% n-update:
ntmp = bsxfun(@times, w,e_lreshaped); % (J x K x I)
z=squeeze(sum(ntmp,2)); % normalization constant (J x I)

t1=e_l.*bsxfun(@minus, alpha', a); %(KxI)
t2=dvec.*log(z); %(JxI)

% approximation of the log-gamma function
t31=gammalogapprox(a)-a.*log(b);
t32=gammalogapprox(alpha)-alpha.*log(beta);
t3=bsxfun(@minus, t31, t32');

l=sum(t1+t3,1)+sum(t2,1); %(1xI)
