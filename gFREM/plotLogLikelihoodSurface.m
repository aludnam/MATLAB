savethis = 0;
p.d=0;
p.lambda = 655; %nm
p.NA = 1.2;
p.pixelsize = 106; %nm
p.sig1=sqrt(2)/2/pi*p.lambda/p.NA/p.pixelsize; %[Zhang 2007]
p.sig2=p.sig1;
p.int1=100;
p.int2=0;
p.bg=0;
x=-25:.01:25;

c1_true = p.d/2; %[0 .5 1 1.5 2 3 4];
c2_true = -c1_true;
q=10;
% c1=[c1_true-q:.2:c1_true+q];
% c2=[c2_true-q:.2:c2_true+q];
c1=[-q:.2:q];
c2=[-q:.2:q];
f1_true = makeGauss(x,c1_true,p.sig1);
f2_true = makeGauss(x,c2_true,p.sig2);
l_true = p.int1*f1_true+p.int2*f2_true+p.bg;
logl=zeros(length(c1),length(c2));

for ii=1:length(c1) 
    f1=makeGauss(x,c1(ii),p.sig1);                % creates PSF (gauss approx)
    for jj=1:length(c2)
        f2=makeGauss(x,c2(jj),p.sig2);            % creates PSF (gauss approx)        
        l=p.int1*f1+p.int2*f2+p.bg;
%         l=p.int1*f1+p.int2*f2;
        logl(ii,jj)=sum(l_true.*log(l));
    end
end
figure; 
% [DX,DY]=gradient(logl, c1(2)-c1(1));
% quiver(c1,c2,DX,DY)
% hold on
% contour(c1,c2,logl, 20)
% colorbar
% set(gca,'dataaspectratio',[1 1 1])
mesh(c1,c2,logl)
hold on
scatter3(c1_true, c2_true, max(logl(:)));
scatter3(c2_true, c1_true, max(logl(:)))
title(['bg=' num2str(p.bg)]);

xlabel('c_1 (center of the source one)')
ylabel('c_2 (center of the source two)')
setfontsizefigure(12)
% colormap(gray)
if savethis 
    name = ['images/LogLikelihoodSurface_d' num2str(d)];
    SaveImageFULL(name)
end