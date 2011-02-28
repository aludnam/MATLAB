function out=bgsubtractbyhand(in,peval)
% out=bgsubtractbyhand(in,peval.bg)
mfprintf(peval.fid,'Background (%g) subtracted and clipped by hand.\n',peval.bg)
out=max(in-peval.bg,eps);