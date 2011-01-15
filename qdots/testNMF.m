
xmat = [15 15 15+p.separ(1) 15;
    15+p.separ(1) 15 15 15;
    15+p.separ(1)/2 15 15+p.separ(1)/2 15;
    15 15 16 15;
    15 15 15+p.separ(1) 16];
for ixmat=1:size(xmat,1)
    x=xmat(ixmat,:);
    wg = makegauss(x, p.s, [p.nx p.ny]);
    w = [reshape(wg, p.nx*p.ny,size(wg,3)), wbg];
    % normalization of all w:
    sumw = sum(w,1);
    w = w./repmat(sumw,n,1); %normalization of each component
    h = h.*repmat(sumw',1,m); %to keep the multiplication equal
    
    x1=repmat(sum(w,1)',1,m);
    %%%%%test
    h=hinit;
    iih=1;
    for ih = 1: maxh
        y1=w'*(v./(w*h));
        h(h_dovec,:)=h(h_dovec,:).*(y1(h_dovec,:))./x1(h_dovec,:);
        h=max(h,eps); % adjust small values to avoid undeflow
        d(k) = ddivergence(v, w*h);
        %%%testing
        dh(ixmat,iih)=d(k);
        iih=iih+1;
        %%%testing
        
        dd(k) = abs(d(k)-d(k-1));
        fprintf('[%g] Ddivergence %g\n',k, d(k))
        k=k+1;
    end
end

plot(dh')
makelegend(xmat,'pos','iter H', 'ddiv');