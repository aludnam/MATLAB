function locerrormeasure=computeerror(xt,yt,xl,yl)
% locerrormeasure=computeerror(xt,yt,xl,yl)
[dist, indexmutual]=locerr(xt, yt, xl, yl);
locerrormeasure=mean(dist);