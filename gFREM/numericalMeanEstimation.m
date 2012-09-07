function I=numericalMeanEstimation(xhires,f1,f2, offset, nphot, osf)
% I=numericalMeanEstimation(xhires,f1,f2, offset, nphot, osf)

if ~exist('offset','var')
    offset=0;
end
if ~exist('nphot','var')
    nphot=1;
end

% individual states of lambda l(:,alpha)=L^alpha_1*q_1+L^alpha_2*q_2
l(:,1)=f1;
l(:,2)=f2;
l(:,3)=(f1+f2);
l(:,4)=zeros(length(f1),1);

if ndims(xhires)==2 %1D vector
    dx=xhires(2)-xhires(1);
else
    dx=xhires(1,2,1)-xhires(1,1,1);
end
N=size(l,2);

precMask=eps;            % Precision for masking significant values in derivative (avoiding division with small numbers) 
prec=eps;                % Precision up to which n_k will be generated                    

for ii=1:N;
    l2d=reshape(l(:,ii),size(xhires,1),size(xhires,2));   
    llowres = nphot*binsumImage(l2d,[osf,osf]);
    lbg(:,ii)=reshape(llowres,numel(llowres),1)+offset;                 % pixelised version and added background in each pixel 
    dl2d=-gradient(reshape(l2d,size(xhires,1),size(xhires,2)),dx);      % Derivative of the intensity of the Poisson 
    dl2dlowres=nphot*binsumImage(dl2d,[osf,osf]);
    dl(:,ii)=reshape(dl2dlowres,numel(dl2dlowres),1);   
end

pcdf=poisscdf(0:50*max(lbg(:)),max(lbg(:)));
nmax = sum(pcdf<1-prec);

n=0:nmax;
lf=size(lbg,1); 
nmat = repmat(n',1,lf);
ln=length(n);


r=zeros(ln,lf,N);
Po=zeros(ln,lf,N);
dlbgMat=zeros(ln,lf,N);

for jj=1:N  % Over 4 states
    lbgMat = repmat(lbg(:,jj)',ln,1);                           % Intensity of the Poisson
    %dlbgMat(:,:,jj)= repmat(gradient(lbg(:,jj),dx)',ln,1);     % <-!!!WRONG!!! direction of derivative     
    dlbgMat(:,:,jj)=repmat(dl(:,jj)',ln,1);

    r(:,:,jj)=(nmat-lbgMat)./lbgMat;                            % Linear ramp with offset and scaling
    Po(:,:,jj)=poissonpdfmulti(n,lbg(:,jj));                    % Poisson distribution
end

rPo=Po.*r;          % Nominator term in the brackets.
sPo=sum(Po,3);      % Denominator term...

It11=sum(rPo(:,:,[1,3]),3).^2;
It22=sum(rPo(:,:,[2,3]),3).^2;
It12=sum(rPo(:,:,[1,3]),3).*sum(rPo(:,:,[2,3]),3);

dl1=dlbgMat(:,:,1);
dl2=dlbgMat(:,:,2);

% Multiplication (masked) with the derivative and division with sPo (to avoid division with small numbers) and sum over n_k and pixels: 

tmp = dl1.^2.*It11;
mask=tmp>precMask;
Ie11=sum(tmp(mask)./sPo(mask),1)';

tmp = dl2.^2.*It22;
mask=tmp>precMask;
Ie22=sum(dl2(mask).^2.*It22(mask)./sPo(mask),1)';

tmp=dl1.*dl2.*It12;
mask=tmp>precMask;
Ie12=sum(dl1(mask).*dl2(mask).*It12(mask)./sPo(mask),1)';

I(1,1)=1/4*Ie11;
I(2,2)=1/4*Ie22;
I(1,2)=1/4*Ie12;
I(2,1)=I(1,2); 