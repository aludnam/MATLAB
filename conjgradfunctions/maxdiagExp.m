function f=maxdiagExp(x,varargin)
% Maximizaton of the diagonality of hte matrix. 
covmat = varargin{1};   % #pix X #pix X #taus 
sizevec = varargin{2};
W = exp(reshape(x, sizevec(1)*sizevec(2), sizevec(3))); % #pix X #comp
% W = max(W,eps);
ntau = size(covmat,3); % #taus (different time delays)
ftmp = zeros(1,ntau);
fcol = zeros(1,sizevec(3));
for tau=1:ntau
    covmattmp = covmat(:,:,tau);
    for col=1:sizevec(3)
        wcol = W(:,col);
        fcol(col) = log(wcol'*covmattmp*wcol);
    end
    ftmp(tau) = 0.5*sum(fcol);
end

f = sum(ftmp);