function [W] = amuse(X)
% BSS using eigenvalue value decomposition
% Program written by A. Cichocki and R. Szupiluk
% 
% X [m x N] matrix of observed (measured) signals,
% W separating matrix,
% y estimated separated sources
% p time delay used in computation of covariance matrices
% optimal time-delay default p= 1
%
% First stage: Standard prewhitening

[m,N]=size(X);
if nargin==1,
 n=m; % 
 end;

Rxx=(X*X')/N;

[Ux,Dx,Vx]=svd(Rxx);
 Dx=diag(Dx);
% n=xxx;
 if n<m, % under assumption of additive white noise and
        %when the number of sources are known or can a priori estimated
  Dx=Dx-real((mean(Dx(n+1:m))));
  Q= diag(real(sqrt(1./Dx(1:n))))*Ux(:,1:n)';
%
else    % under assumption of no additive noise and when the 
        % number of sources is unknown
   n=max(find(Dx>1e-199)); %Detection the number of sources
   Q= diag(real(sqrt(1./Dx(1:n))))*Ux(:,1:n)';
end;
%
% else    %assumes no noise
%    Q=inv(sqrtm(Rxx));
% end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second stage: Fast separation using sorting EVD
% notation the same as used in the Chapter 4
Xb=Q*X;
p=1;    
% paramter p can take here value different than 1
% for example -1 or 2.
N=max(size(Xb));
Xb=Xb-kron(mean(Xb')',ones(1,N));

Rxbxbp=(Xb(:,1:N-1)*Xb(:,2:N)')/(N-1);
Rxbxbp= Rxbxbp+Rxbxbp';
[Vxb Dxb]=eig(Rxbxbp);
[D1 perm]=sort(diag(Dxb));
D1=flipud(D1);
Vxb=Vxb(:,flipud(perm));
W = Vxb'*Q;
%y = Vxb' * x1;


