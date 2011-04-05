function [meanf,varf,entf,kurt]=computeMeanVarFromPermutations(f,k)
% [meanf,varf]=computeMeanVarFromPermutations(f,k)
% Computes mean and variance of the all combination of k columns f (f: txn)
n=size(f,2);
choosevec=nchoosek(1:n, k); % (1:nchoosek(n,k)xk) indices to be taken (all combinations) 
nk=nchoosek(n,k);
meanf=zeros(1,nk);
varf=zeros(1,nk);
entf=zeros(1,nk);
kurt=zeros(1,nk);
for ii=1:nk
    ftmp = sum(f(:,choosevec(ii,:)),2); %tx1
    meanf(ii)=mean(ftmp(:));
    varf(ii)=var(ftmp(:));
%     ftmppos=ftmp+min(ftmp); %to make it positiove everywhere
%     [ftmphist,X] = hist(ftmppos,100);
%     entf(ii)=computeEntropy(ftmphist);
    entf(ii)=histogrametropy(ftmp);
    kurt(ii)=kurtosis(ftmp)-3;
end