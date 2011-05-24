showim=1;
% positions
l1=0;
% l2=[.1:.1:2];
l2=[0:.5:4];
% sigmas
sig1=1;
sig2=2;
% intensities
int1= 1;
int2=10;

x=[-5:.01:20];

for ii=1:length(l2)
    f1=makeGauss(x,l1,sig1);
    f2=makeGauss(x,l2(ii),sig2);
    
    y(ii)=gFREMfunction(x,f1,f2, int1, int2,showim);
end

y1=gFREMfunction(x,f1,f2,int1,0);
y2=gFREMfunction(x,f1,f2,0,int2);