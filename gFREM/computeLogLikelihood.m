function [logl, c1, c2, p] = computeLogLikelihood(d,bg)
% [logl, c1, c2, p] = computeLogLikelihood(d,bg)
p.d=d;
p.bg=bg;
p.lambda = 655; %nm
p.NA = 1.2;
p.pixelsize = 106; %nm
p.sig1=sqrt(2)/2/pi*p.lambda/p.NA/p.pixelsize; %[Zhang 2007]
p.sig2=p.sig1;
p.int1=750;
p.int2=750;

x=-25:.1:25;

p.c1_true = p.d/2; %[0 .5 1 1.5 2 3 4];
p.c2_true = -p.c1_true;
q=8;
qstep = .5;
% c1=[p.c1_true-q:.2:p.c1_true+q];
% c2=[p.c2_true-q:.2:p.c2_true+q];
c1=[-q:qstep:q];
c2=[-q:qstep:q];
f1_true = makeGauss(x,p.c1_true,p.sig1);
f2_true = makeGauss(x,p.c2_true,p.sig2);
l_true = p.int1*f1_true+p.int2*f2_true+p.bg;
logl=zeros(length(c1),length(c2));

for ii=1:length(c1) 
    f1=makeGauss(x,c1(ii),p.sig1);                % creates PSF (gauss approx)
    for jj=1:length(c2)
        f2=makeGauss(x,c2(jj),p.sig2);            % creates PSF (gauss approx)        
        l=p.int1*f1+p.int2*f2+p.bg;       
        logl(ii,jj)=sum(l_true.*log(l)-l-factorialapprox(l_true));
    end
end
