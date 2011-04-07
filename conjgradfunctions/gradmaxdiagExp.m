function gf=gradmaxdiagExp(x,varargin)
% Maximizaton of the diagonality of hte matrix.
covmat = varargin{1};   % #pix X #pix X #taus
sizevec = varargin{2};
W = exp(reshape(x, sizevec(1)*sizevec(2), sizevec(3)));
ntau = size(covmat,3); % #taus (different time delays)
gftmp = zeros(ntau,prod(sizevec));
fcol = zeros(1,sizevec(3));
for tau=1:ntau
    covmattmp = covmat(:,:,tau);
    for col=1:sizevec(3)
        wcol = W(:,col);
        fcol(col) = (wcol'*covmattmp*wcol);
    end
    Q = diag(1./fcol);
    gftmp(tau,:) = reshape((Q*W'*0.5*(covmattmp+covmattmp'))',1,prod(sizevec));        
end
gf = sum(gftmp,1).*exp(x);