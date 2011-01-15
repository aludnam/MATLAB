function [Gn, pi1,pi2, am_pix] = icacalcpi( W, Q, A );
% Performance indices for extraction/separation
% (SCruces 01Jan03)
%
% These performance indices resembles Amari's index 
% and work properly when the sources are normalized 
% to unit variance.
if isempty(W) | (size(W,2) ~= size(A,1))
   disp('First press button "RUN ALGORITHM"');
   Gn = []; pi1 = []; pi2 = [];
   return 
end

G = ( W * Q * A );
[Nrows, Ncols] = size( G );

% pi1: Gives an idea of the capability of extraction that the algorithm
% achieved. Computes some normalized mean interference due to the other sources
% that remains in the extracted/separated signals. This index ideally 
% should be zero. However, it does not take into account if the 
% extracted sources are different or not (this can be controled with the 
% performance index pi2).


D=pinv(diag(max(abs(G.'))));
Gn =D*abs(G);
pi1=(sum(sum(Gn))-Nrows)/(Nrows*Ncols-Nrows);

% pi2: Indicates the divergence from Unitarity/Semi-Unitarity 
% after the normalization of the rows of G to unit norm to remove
% any possible scaling indeterminacy. Idealy should be zero.

GG=abs(G*G');
D=pinv(diag(sqrt(diag(GG))));
GGn=D*GG*D;
if Nrows>1
    pi2=(sum(sum(GGn))-Nrows)/(Nrows*(Nrows-1));
else
    pi2=0;
end

% Amari performance index
am_pix = sum(sum(abs(G)./repmat(max(abs(G + eps),[],2),1,Ncols) + abs(G)./repmat(max(abs(G + eps),[],1),Nrows,1)))/(2*Ncols) - 1;

