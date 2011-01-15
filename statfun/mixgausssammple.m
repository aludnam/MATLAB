function x=mixgausssammple(a,b,p,N,M)
% x=mixgausssammple(a,b,p,N,M)
% N x M sample from mixture of two gamma distribution
% a,b: shape and scale parameters of the gamma
% p: mixing coeficient of the mixture

n1=binornd(N,p,M,1);
for ii=1:M
    x1=gamrnd(a(1),b(1),1,n1(ii));
    x2=gamrnd(a(2),b(2),1,N-n1(ii));
    xtmp(ii,:)=[x1,x2];
end
indexperm=randperm(N);
x=xtmp(:,indexperm); %random permutaiton of the xtmp
