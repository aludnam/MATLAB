function f=maxdiagnoExp(x,varargin)
% Maximizaton of the diagonality of hte matrix. 
covmat = varargin{1};   % #pix X #pix X #taus 
sizevec = [size(covmat,1), length(x)/size(covmat,1)];
W = pinv(exp(reshape(x, sizevec(1), sizevec(2))))'; % #pix X #comp
% W = max(W,eps);
ntau = size(covmat,3); % #taus (different time delays)
ftmp = zeros(1,ntau);
for tau=1:ntau
    covmattmp = covmat(:,:,tau);
    M = W'*covmattmp*W;
%     ftmp(tau) = 0.5*(sum(log(diag(M)))-log(abs(det(M))));
ftmp(tau) = 0.5*abs(sum(log(abs(diag(M))))-log(abs(det(M))));
end
f = sum(ftmp);