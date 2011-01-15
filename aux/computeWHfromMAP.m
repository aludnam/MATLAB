function [Wxk,Hkt,centers,Vxtpix, Vxtpixbg] = computeWHfromMAP(res,peval)
% [Wxk,Hkt,centers,Vxtpix, Vxtpixbg] = computeWHfromMAP(res,peval)

[Wxk,Hkt,centers,Vxtpix]=reshapeGaP(res.hvec,res.cxcy,peval);
Vxtpixbg=reshape(Wxk*Hkt,peval.nx,peval.ny,peval.nt)+peval.bg;