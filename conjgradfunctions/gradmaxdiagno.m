function gf=gradmaxdiagno(x,varargin)
% Maximizaton of the diagonality of hte matrix. 
covmat = varargin{1};   % #pix X #pix X #taus 
sizevec = varargin{2};
W = reshape(x, sizevec(1)*sizevec(2), sizevec(3));
ntau = size(covmat,3); % #taus (different time delays)
gftmp = zeros(ntau,prod(sizevec));
for tau=1:ntau
    covmattmp = covmat(:,:,tau);
    P = W'*covmattmp*W;
    Q = diag(1./diag(P));
    R = inv(P);
    gfr = (Q'*W' - (W*R)')*0.5*(covmattmp+covmattmp');
    gftmp(tau,:) = reshape(gfr',1,prod(sizevec));
end

gf = sum(gftmp,1);