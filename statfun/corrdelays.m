function C=corrdelays(A, tdelay)
% C=corrdelays(A, tdelays)
% Computes correlation matrices of MxN matrix A (rows are time series) for
% time delayes specified in the 1xT vector TDELAY

[m,n]=size(A);
t=length(tdelay);
C=zeros(m,m,t);
for indext = 1:t
    Ashift = (circshift(A',tdelay(indext)))';
    C(:,:,indext)=corr(A',Ashift');    
end