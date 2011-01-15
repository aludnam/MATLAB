function [Wxk,Hkt,centers,Wxkpix]=reshapeGaP(x1,x2,peval)
% [Wxk,Hkt,centers,Wxkpix]=reshapeGaP(x1,x2,peval);
Hkt=exp(reshape(x1, peval.ncomp, peval.nt));
%centers=reshape(x2, peval.ncomp,2)+1; %plus one because of the dip_image notation
centers=reshape(x2, peval.ncomp,2);
normconst = 1./(peval.sigmaPSF.^2*2*pi);
Wxkpix=gauss2dmultislice([peval.nx, peval.ny, peval.ncomp], centers+1, peval.sigmaPSF, normconst);
Wxk=reshape(Wxkpix,peval.nx*peval.ny, peval.ncomp);

%just to make sure W is normalized 
Wxk=Wxk./repmat(sum(Wxk,1), peval.nx*peval.ny,1);
Wxkpix=reshape(Wxk,peval.nx,peval.ny,peval.ncomp);

