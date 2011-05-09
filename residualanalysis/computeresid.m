function [resid, resid_norm, Wxk,Hkt,centers,Vxtpix, Vxtpixbg] = computeresid(res,dpixc,peval)
% [resid, resid_norm, Wxk,Hkt,centers,Vxtpix, Vxtpixbg] =
% computeresid(res,dpixc,peval)
[Wxk,Hkt,centers,Vxtpix]=reshapeGaP(res.hvec,res.cxcy,peval);
Vxtpixbg=reshape(Wxk*Hkt,peval.nx,peval.ny,peval.nt)+peval.bg;
resid=(Vxtpixbg-dpixc);
resid_norm=resid./sqrt(Vxtpixbg);
