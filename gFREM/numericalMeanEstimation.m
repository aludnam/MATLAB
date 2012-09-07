function I=numericalMeanEstimation(xhires,f1,f2, offset, nphot, osf)
% I=numericalMeanEstimation(x,f1,f2, offset, nphot)

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

%precMask=10^-8;             % Precision for masking intenisty images     
precMask = 10^-8;%eps; 
prec=10^-10;                % Precision up to which n_k will be generated                    

for ii=1:N;
    l2d=reshape(l(:,ii),size(xhires,1),size(xhires,2));   
    llowres = nphot*(binsumImage(l2d,[osf,osf]));
    lbg(:,ii)=reshape(llowres,numel(llowres),1)+offset; % pixelised version and added background in each pixel 
    dl2d=-gradient(reshape(l2d,size(xhires,1),size(xhires,2)),dx);     % Derivative of the intensity of the Poisson 
    dl2dlowres=binsumImage(dl2d,[osf,osf]); %
    dl(:,ii)=reshape(dl2dlowres,numel(dl2dlowres),1);
    %dlbgMat(:,:,ii)=repmat(reshape(dl2dlowres,1,numel(dl2dlowres)),ln,1);
end

%lbg=l+offset;               % adding background
pcdf=poisscdf(0:50*max(lbg(:)),max(lbg(:)));
nmax = sum(pcdf<1-prec);

n=0:nmax;
lf=size(lbg,1); 
nmat = repmat(n',1,lf);
ln=length(n);


r=zeros(ln,lf,N);
Po=zeros(ln,lf,N);
dlbgMat=zeros(ln,lf,N);
maskdl=zeros(ln,lf,N);
for jj=1:N  % Over 4 states
    lbgMat = repmat(lbg(:,jj)',ln,1);                           % Intensity of the Poisson
    maskdl(:,:,jj)=lbgMat>precMask;                             % Restriction to non-zero lambda (over space) to avoid division with small number
    %dlbgMat(:,:,jj)= repmat(gradient(lbg(:,jj),dx)',ln,1);     % <-!!!WRONG!!! direction of derivative     
    dlbgMat(:,:,jj)=repmat(dl(:,jj)',ln,1);

    r(:,:,jj)=(nmat-lbgMat)./lbgMat;                            % Linear ramp with offset and scaling
    Po(:,:,jj)=poissonpdfmulti(n,lbg(:,jj));                    % Poisson distribution
end

rPo=Po.*r.*maskdl;
sPo=sum(Po,3);      % Denominator term...
% sPo=sum(Po(:,:,1:3),3);

% It11=dlbgMat(:,:,1).^2.*sum(rPo(:,:,[1,3]),3).^2;
% It22=dlbgMat(:,:,2).^2.*sum(rPo(:,:,[2,3]),3).^2;
% It12=dlbgMat(:,:,1).*dlbgMat(:,:,2).*sum(rPo(:,:,[1,3]),3).*sum(rPo(:,:,[2,3]),3);

It11=sum(rPo(:,:,[1,3]),3).^2;
It22=sum(rPo(:,:,[2,3]),3).^2;
It12=sum(rPo(:,:,[1,3]),3).*sum(rPo(:,:,[2,3]),3);
% summation over n_k:
mask = abs(sPo)<precMask; % to avoid division of small numberss...
sPo(mask)=1;It11(min(mask,It11<precMask))=0;It22(min(mask,It22<precMask))=0;It12(min(mask,It12<precMask))=0;

Ie11=sum(It11./sPo,1)';
Ie22=sum(It22./sPo,1)';
Ie12=sum(It12./sPo,1)';
% multiplication with the derivative and sum over pixels: 
mdl1=abs(dl(:,1))>precMask;
mdl2=abs(dl(:,2))>precMask;
mdl3=min(mdl1,mdl2);
I(1,1)=1/4*sum(dl(mdl1,1).^2.*Ie11(mdl1));
I(2,2)=1/4*sum(dl(mdl2,2).^2.*Ie22(mdl2));
I(1,2)=1/4*sum(dl(mdl3,1).*dl(mdl3,2).*Ie12(mdl3));
I(2,1)=I(1,2); 

% maskdl=abs(dl)>eps;
% I(1,1)=1/4*sum(dl(maskdl(:,1),1).^2.*Ie11(maskdl(:,1)')');
% I(2,2)=1/4*sum(dl(maskdl(:,2),2).^2.*Ie22(maskdl(:,2)')');
% mdlcom=min(maskdl(:,1),maskdl(:,2));
% I(1,2)=1/4*sum(dl(mdlcom,1).*dl(mdlcom,2).*Ie12(mdlcom)');
% I(2,1)=I(1,2); 

% dxlowres=dx*osf;
% mask = abs(It11)>precMask;
% I(1,1)=1/4*dxlowres*trapz(It11(mask)./sPo(mask));
% mask = abs(It22)>precMask;
% I(2,2)=1/4*dxlowres*trapz(It22(mask)./sPo(mask));
% mask = abs(It12)>precMask;
% I(1,2)=1/4*dxlowres*trapz(It12(mask)./sPo(mask));
% I(2,1)=I(1,2); 