function [intf, intx]=pixelize1d(x,f)
% [intf, intx]=pixelize1d(x,f)
% Integrates over integer values in x (pixels) the funciton f. 

xfl=floor(x);
indeces=[0,(find(diff(xfl))),length(x)];

for ii=1:(xfl(end)+abs(xfl(1)))
    indbeg=indeces(ii)+1;
    indend=indeces(ii+1);
    
    intf(ii)=trapz(x(indbeg:indend),f(indbeg:indend));
    intx(ii)=xfl(indbeg)+.5;
end

