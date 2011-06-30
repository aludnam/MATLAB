% This computes FREM for different number of states the sources can get
% into (uniformly distributed over these states).
savethis =0;

int1_multi{1}=[500];
% int1_multi{2}=[500];
% int1_multi{3}=[1000];
% int1_multi{4}=[5000];

int2_multi = int1_multi;



% int1_multi{2}=[0 1];
% int2_multi{2}=[0 1];
% int1_multi{3}=[0 .5 1];
% int2_multi{3}=[0 .5 1];
% int1_multi{4}=[0:.1:1];
% int2_multi{4}=[0:.1:1];

tau1_vec = [0 .2 .5 .9];

tau2_vec = tau1_vec;
clear int1_mat;
int1_mat{1}=[0 1];
clear('y_multi')

% positions of sources
l1=0;
l2=0:.2:10;
x=-7:.01:16;
% x=-10:.01:20;

p.offset = 100;

% sigma
p.lambda = 655; %nm
p.NA = 1.2;
p.pixelsize = 106; %nm
p.sig1=sqrt(2)/2/pi*p.lambda/p.NA/p.pixelsize; %[Zhang 2007]
p.sig2=p.sig1;

probfunction='binomial';
nameappendix = 'Ind';
correctIntensity = 0;
if correctIntensity
    nameappendix = [nameappendix 'CorrInt'];
end
pixelizeversion = 0;
vard=zeros(length(l2), length(int1_multi));

for mm=1:length(int1_multi)   
    tau= cat(1,tau1_vec(mm),tau2_vec(mm));
    int_vec=cat(1,int1_multi{mm},int2_multi{mm});
    [pint, int_out]=generateDistribution(int_vec,tau, probfunction, correctIntensity);
    [vard(:,mm), I3d(:,:,:,mm)]=computeSeparationVariance(x,l1,l2,[p.sig1,p.sig2],int_vec, pint, pixelizeversion, p.offset);    
%     if length(int1_multi{mm})==1
%         % Integrating out
%         [vardintout(:,mm), Iintout(:,:,:,mm)]=computeSeparationVarianceIntOut(x,l1,l2,[p.sig1,p.sig2],int_vec, pixelizeversion, p.offset);
%     end
end

% plotledacos

