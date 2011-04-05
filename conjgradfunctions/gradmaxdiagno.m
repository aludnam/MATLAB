function gf=gradmaxdiagno(x,varargin)
% Maximizaton of the diagonality of hte matrix. 
covmat = varargin{1};   % #pix X #pix X #taus 
sizevec = varargin{2};
W = reshape(x, sizevec(1)*sizevec(2), sizevec(3));
ntau = size(covmat,3); % #taus (different time delays)
gftmp = zeros(ntau,prod(sizevec));
for tau=1:ntau
    covmattmp = covmat(:,:,tau);
    Q = diag(1./diag(W'*covmattmp*W));
    gftmp(tau,:) = reshape(Q'*W'*0.5*(covmattmp+covmattmp'),1,prod(sizevec));
end

gf = sum(gftmp);