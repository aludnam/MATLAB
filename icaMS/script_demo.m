% Demo script
% version 1.1


% Generate data
K_true=6;        % Number of sources
M=10;            % Number of observations
N=1000;          % Number of samples
interval=[1:M];  % Number of IC components interval to search

disp(sprintf('True number of sources: %d\n',K_true));

% Generate mixed signals from artificial super Gaussian sources and mixing matrix
Snorm=randn(N,K_true);
Atrue=randn(M,K_true);
Strue=filter([1 .4], 1,sinh(Snorm)/.1)';
X=Atrue*Strue;

% remove mean
X=X-repmat(mean(X,2),1,N);




% - 1th Part - Finding best model wrt. number of components using BIC

% BIC
[P,logP]=icaMS_bic(X,interval,1);

[most_prop,most_prop_K]=max(P);
most_prop_K=interval(most_prop_K);

disp(sprintf('\nEstimated components %d with P=%0.2f\n',most_prop_K,most_prop));


% Draw
figure(1)
clf
bar(interval,P);
drawnow;




% - 2nd Part - finding ICA components using BIC estimat of K (number of source components)

disp(sprintf('PCA and MS %d components:',most_prop_K));

% PCA
[U,D,V]=svd(X',0);

Xica=( U(:,1:most_prop_K)*D(1:most_prop_K,1:most_prop_K) )';


% ICA
[S,Aica,ll]=icaMS(Xica);
A=V(:,1:most_prop_K)*Aica;



