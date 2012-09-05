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
precMask = eps; 
prec=10^-10;                % Precision up to which n_k will be generated                    

for ii=1:size(l,2);
    l2d=reshape(l(:,ii),size(xhires,1),size(xhires,2));   
    llowres = nphot*normalize(imresize(l2d,1/osf));
    lbg(:,ii)=reshape(llowres,numel(llowres),1)+offset; % pixelised version and added background in each pixel 
    dl2d=gradient(reshape(l2d,size(xhires,1),size(xhires,2)),dx);     % Derivative of the intensity of the Poisson 
    dl2dlowres=osf^2*imresize(dl2d,1/osf); % the factor osf^2 is there for keeping the total intensity (integrated over the whole image) unchanged
    dl(:,ii)=reshape(dl2dlowres,numel(dl2dlowres),1);
    %dlbgMat(:,:,ii)=repmat(reshape(dl2dlowres,1,numel(dl2dlowres)),ln,1);
end

%lbg=l+offset;               % adding background
pcdf=poisscdf(0:50*max(lbg(:)),max(lbg(:)));
nmax = sum(pcdf<1-prec);

n=0:nmax;
ln=length(n);
lf=size(lbg,1);
nmat = repmat(n',1,lf);
r=zeros(ln,lf,N);
Po=zeros(ln,lf,N);
dlbgMat=zeros(ln,lf,N);
maskdl=zeros(ln,lf,N);
for jj=1:N  % Over 4 states
    lbgMat = repmat(lbg(:,jj)',ln,1);                           % Intensity of the Poisson
    maskdl(:,:,jj)=lbgMat>precMask;                             % Restriction to non-zero lambda (over space)
    %dlbgMat(:,:,jj)= repmat(gradient(lbg(:,jj),dx)',ln,1);     % <-!!!WRONG!!! direction of derivative     
    dlbgMat(:,:,jj)=repmat(dl(:,jj)',ln,1);
    r(:,:,jj)=(nmat-lbgMat)./lbgMat;                            % Linear ramp with offset and scaling
    Po(:,:,jj)=poissonpdfmulti(n,lbg(:,jj));                    % Poisson distribution
end

rPo=Po.*r.*maskdl;
sPo=sum(Po,3);      % Denominator term...
% sPo=sum(Po(:,:,1:3),3);


It11=dlbgMat(:,:,1).^2.*sum(rPo(:,:,[1,3]),3).^2;
It22=dlbgMat(:,:,2).^2.*sum(rPo(:,:,[2,3]),3).^2;
It12=dlbgMat(:,:,1).*dlbgMat(:,:,2).*sum(rPo(:,:,[1,3]),3).*sum(rPo(:,:,[2,3]),3);

dxlowres=dx*osf;
mask = abs(It11)>precMask;
I(1,1)=1/4*dxlowres*trapz(It11(mask)./sPo(mask));
mask = abs(It22)>precMask;
I(2,2)=1/4*dxlowres*trapz(It22(mask)./sPo(mask));
mask = abs(It12)>precMask;
I(1,2)=1/4*dxlowres*trapz(It12(mask)./sPo(mask));
I(2,1)=I(1,2); 