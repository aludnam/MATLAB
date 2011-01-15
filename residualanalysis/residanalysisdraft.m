Wxkpixbg=reshape(Wxk*Hkt,peval.nx,peval.ny,peval.nt)+peval.bg;
resid=(Wxkpixbg-dpixc);
resid_norm=resid./Wxkpixbg;
rnm=mean(resid_norm,3);