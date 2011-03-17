function [meanf,varf,entf]=computeMeanVarFromPermutations(f,k)
% [meanf,varf]=computeMeanVarFromPermutations(f,k)
% Computes mean and variance of the all combination of k columns f (f: txn)
n=size(f,2);
choosevec=nchoosek(1:n, k); % (1:nchoosek(n,k)xk) indices to be taken (all combinations) 
for ii=1:nchoosek(n,k)
    ftmp = sum(f(:,choosevec(ii,:)),2); %tx1
    meanf(ii)=mean(ftmp(:));
    varf(ii)=var(ftmp(:));
%     ftmppos=ftmp+min(ftmp); %to make it positiove everywhere
%     [ftmphist,X] = hist(ftmppos,100);
%     entf(ii)=computeEntropy(ftmphist);
    entf(ii)=histogrametropy(ftmp);
end