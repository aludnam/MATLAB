k=2;
nx=size(Vxt,1);
nt=size(Vxt,2);
Hkte=rand(k,nt);
Wxke=rand(nx,k);
sumW=sum(exp(Wxke),1);
Wxken=Wxke+log(1./repmat(sumW,nx,1)); %to make sum(exp(Wxken))=ones(1,k)
Hkte= Hkte-log(1./repmat(sumW',1,nt)); %to keep exp(Wxken)*exp(Hkte) constant
%background:
[out, bg, bg_im]=backgroundoffset(Vxt, 'no', 5, 30, 8);
Wxk_fix=1/nx*ones(nx,1);    
Hkt_fix=nx*bg*ones(1,nt);

options = zeros(1,18);
options (1)=1; %to display error values
options (7)=1;
options(9)=0; %to check gradient
options(14)=1; %maximum number of iterations

for ii=1:200
    [rHkte, options, flog, pointlog] = conjgrad('ddivHexp', reshape(Hkte,1,nt*k), options, 'gradddivHexp', Vxt, Wxken, Wxk_fix, Hkt_fix);
    Hkte=reshape(rHkte,k,nt);
    dH(ii)=ddivHexp(reshape(Hkte,1,nt*k),Vxt, Wxken, Wxk_fix, Hkt_fix);
    [rWxke, options, flog, pointlog] = conjgrad('ddivWexp', reshape(Wxken,1,nx*k), options, 'gradddivWexp', Vxt, Hkte, Wxk_fix ,Hkt_fix);
    Wxke=reshape(rWxke,nx,k);
    %normalization of Wxk:
    sumW=sum(exp(Wxke),1);
    Wxken=Wxke+log(1./repmat(sumW,nx,1)); %to make sum(exp(Wxken))=ones(1,k)
    Hkte= Hkte-log(1./repmat(sumW',1,nt)); %to keep exp(Wxken)*exp(Hkte) constant
    dW(ii)=ddivWexp(Wxken,Vxt, Hkte, Wxk_fix, Hkt_fix);
    
end