%bd = fitH (x,m,n)
%computes  beta distrib in points in x for parameters m and n
function bd = fitH(x,p)
nf = gamma(p(1))*gamma(p(2))/gamma(p(2)+p(1));
bd = 1/nf * x.^(p(1)-1).*(1-x).^(p(2)-1); %gamma distrib...
