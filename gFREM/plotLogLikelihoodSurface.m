sig1=1;
sig2=1;
int1=100;
int2=100;
x=-10:.01:25;
c1=-5:.1:5;
c2=-5:.1:5;
logl=zeros(length(c1),length(c2));
for ii=1:length(c1) 
    f1=makeGauss(x,c1(ii),sig1);            % creates PSF (gauss approx)
    for jj=1:length(c2)
        f2=makeGauss(x,c2(jj),sig2);            % creates PSF (gauss approx)        
        l=int1*f1+int1*f2;
        logl(ii,jj)=sum(l.*log(l)-l-factorialapprox(l));
    end
end


ifcontourf(logl)
colorbar
set(gca,'dataaspectratio',[1 1 1])
xlabel('c_1 (center of the source one)')
ylabel('c_2 (center of the source two)')
setfontsizefigure(12)
colormap(gray)
if savethis 
    SaveImageFULL('images/LogLikelihoodSurface')
end